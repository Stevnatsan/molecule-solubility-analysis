(function () {
  "use strict";

  const $ = (id) => document.getElementById(id);
  const NS = "http://www.w3.org/2000/svg";
  const minus = (s) => String(s).replace(/^-/, "−");
  const fmt = (v, digits = 2) => minus(Number(v).toFixed(digits));
  const sig = (v) => Number(v.toPrecision(2)).toLocaleString("en-US", { maximumFractionDigits: 3 });

  function el(tag, attrs = {}, text) {
    const node = tag.startsWith("svg:") ? document.createElementNS(NS, tag.slice(4)) : document.createElement(tag);
    for (const [k, v] of Object.entries(attrs)) node.setAttribute(k, v);
    if (text != null) node.textContent = text;
    return node;
  }

  const data = JSON.parse($("dashboard-data").textContent);
  const structures = JSON.parse($("structure-data").textContent);
  const cols = data.columns;
  const mols = data.molecules.map((row, id) => {
    const m = { id };
    cols.forEach((c, i) => { m[c] = row[i]; });
    m.residual = m.logS - m.pred;
    m.bad = Math.abs(m.residual) > 1;
    return m;
  });
  const overallRmse = data.metrics.randomForest.rmse;
  const familyByName = Object.fromEntries(data.families.map((f) => [f.name, f]));

  const AXES = {
    logP: { button: "logP", title: "logP (Crippen) — higher = oilier", get: (m) => m.logP },
    MW: { button: "Weight", title: "Molecular weight (g/mol)", get: (m) => m.MW },
    TPSA: { button: "Polarity", title: "Polar surface area (Å²)", get: (m) => m.TPSA },
    pred: { button: "Predicted", title: "Model's predicted logS", get: (m) => m.pred },
  };

  const state = { x: "logP", family: null, selected: null };

  /* ---------- KPIs ---------- */
  const organo = familyByName["Organophosphate"];
  [
    ["Molecules analysed", data.metrics.n.toLocaleString("en-US"), "", ""],
    ["Model error (cross-validated)", fmt(overallRmse), "log units", "kpi--accent"],
    ["Published ESOL equation", fmt(data.metrics.esolEquation.rmse), "log units", ""],
    ["Organophosphate error", fmt(organo.rmse), `${(organo.rmse / overallRmse).toFixed(1)}× average`, "kpi--warn"],
  ].forEach(([label, value, unit, cls]) => {
    const box = el("div", { class: `kpi ${cls}`.trim() });
    const dd = el("dd", {}, value);
    if (unit) dd.append(el("small", {}, unit));
    box.append(el("dt", {}, label), dd);
    $("kpis").append(box);
  });

  /* ---------- controls ---------- */
  for (const [key, axis] of Object.entries(AXES)) {
    const b = el("button", { type: "button", "aria-pressed": String(key === state.x) }, axis.button);
    b.addEventListener("click", () => { state.x = key; syncControls(); drawScatter(); });
    b.dataset.axis = key;
    $("x-axis").append(b);
  }
  const allChip = el("button", { type: "button", class: "chip", "aria-pressed": "true" }, "All");
  allChip.addEventListener("click", () => setFamily(null));
  $("families").append(allChip);
  for (const f of data.families) {
    const chip = el("button", { type: "button", class: "chip", "aria-pressed": "false" }, `${f.name} (${f.n})`);
    chip.dataset.family = String(f.id);
    chip.addEventListener("click", () => setFamily(state.family === f.id ? null : f.id));
    $("families").append(chip);
  }
  function setFamily(id) { state.family = id; syncControls(); drawScatter(); }
  function syncControls() {
    for (const b of $("x-axis").children) b.setAttribute("aria-pressed", String(b.dataset.axis === state.x));
    for (const c of $("families").children) {
      const on = c === allChip ? state.family === null : Number(c.dataset.family) === state.family;
      c.setAttribute("aria-pressed", String(on));
    }
    for (const r of $("family-bars").querySelectorAll("button")) {
      r.setAttribute("aria-pressed", String(Number(r.dataset.family) === state.family));
    }
  }

  const names = $("molecule-names");
  for (const m of [...mols].sort((a, b) => a.name.localeCompare(b.name))) names.append(el("option", { value: m.name }));
  function runSearch() {
    const q = $("search").value.trim().toLowerCase();
    if (!q) { $("search-hint").textContent = ""; return; }
    const hit = mols.find((m) => m.name.toLowerCase() === q) || mols.find((m) => m.name.toLowerCase().includes(q));
    if (!hit) { $("search-hint").textContent = `No molecule matching “${$("search").value.trim()}” in this dataset.`; return; }
    $("search-hint").textContent = "";
    if (state.family !== null && !hit.families.includes(state.family)) setFamily(null);
    select(hit.id);
  }
  $("search").addEventListener("change", runSearch);
  $("search").addEventListener("keydown", (e) => { if (e.key === "Enter") { e.preventDefault(); runSearch(); } });

  /* ---------- scatter ---------- */
  const svg = $("scatter");
  const frame = $("scatter-frame");
  const tooltip = $("tooltip");
  let plotted = []; // points that respond to hover/click (only the highlighted family when one is selected)
  const coords = new Map(); // id -> pixel position for every molecule

  function niceTicks(min, max, count) {
    const raw = (max - min) / count;
    const mag = 10 ** Math.floor(Math.log10(raw));
    const r = raw / mag;
    const step = (r >= 7.5 ? 10 : r >= 3.5 ? 5 : r >= 1.5 ? 2 : 1) * mag;
    const out = [];
    for (let v = Math.ceil(min / step) * step; v <= max + 1e-9; v += step) out.push(Number(v.toFixed(6)));
    return out;
  }

  function fitLine(xs, ys) {
    const n = xs.length, mx = xs.reduce((a, b) => a + b) / n, my = ys.reduce((a, b) => a + b) / n;
    let sxy = 0, sxx = 0, syy = 0;
    for (let i = 0; i < n; i++) { sxy += (xs[i] - mx) * (ys[i] - my); sxx += (xs[i] - mx) ** 2; syy += (ys[i] - my) ** 2; }
    return { slope: sxy / sxx, intercept: my - (sxy / sxx) * mx, r: sxy / Math.sqrt(sxx * syy) };
  }

  const svgSize = () => { const r = svg.getBoundingClientRect(); return { w: r.width, h: r.height }; };

  function drawScatter() {
    const { w, h } = svgSize();
    const pad = { l: 46, r: 14, t: 12, b: 42 };
    const axis = AXES[state.x];
    const xs = mols.map(axis.get), ys = mols.map((m) => m.logS);
    let x0 = Math.min(...xs), x1 = Math.max(...xs);
    let y0 = Math.min(...ys), y1 = Math.max(...ys);
    if (state.x === "pred") { x0 = y0 = Math.min(x0, y0); x1 = y1 = Math.max(x1, y1); }
    const xp = (x1 - x0) * 0.04, yp = (y1 - y0) * 0.05;
    x0 -= xp; x1 += xp; y0 -= yp; y1 += yp;
    const sx = (v) => pad.l + ((v - x0) / (x1 - x0)) * (w - pad.l - pad.r);
    const sy = (v) => h - pad.b - ((v - y0) / (y1 - y0)) * (h - pad.t - pad.b);

    svg.replaceChildren();
    svg.setAttribute("viewBox", `0 0 ${w} ${h}`);
    const axes = el("svg:g", { class: "axis" });
    for (const t of niceTicks(x0, x1, Math.max(3, Math.round(w / 90)))) {
      axes.append(el("svg:line", { x1: sx(t), x2: sx(t), y1: pad.t, y2: h - pad.b }));
      axes.append(el("svg:text", { x: sx(t), y: h - pad.b + 16, "text-anchor": "middle" }, minus(t)));
    }
    for (const t of niceTicks(y0, y1, 6)) {
      axes.append(el("svg:line", { x1: pad.l, x2: w - pad.r, y1: sy(t), y2: sy(t) }));
      axes.append(el("svg:text", { x: pad.l - 8, y: sy(t) + 4, "text-anchor": "end" }, minus(t)));
    }
    axes.append(el("svg:text", { class: "axis-title", x: (pad.l + w - pad.r) / 2, y: h - 6, "text-anchor": "middle" }, axis.title));
    axes.append(el("svg:text", { class: "axis-title", x: -(pad.t + h - pad.b) / 2, y: 13, transform: "rotate(-90)", "text-anchor": "middle" }, "Measured logS (mol/L, log₁₀)"));
    svg.append(axes);

    const clip = el("svg:clipPath", { id: "plot-area" });
    clip.append(el("svg:rect", { x: pad.l, y: pad.t, width: w - pad.l - pad.r, height: h - pad.t - pad.b }));
    svg.append(clip);
    const fit = fitLine(xs, ys);
    if (state.x === "pred") {
      svg.append(el("svg:line", { class: "guide", "clip-path": "url(#plot-area)", x1: sx(x0), y1: sy(x0), x2: sx(x1), y2: sy(x1) }));
    } else if (state.x === "logP") {
      svg.append(el("svg:line", { class: "guide", "clip-path": "url(#plot-area)", x1: sx(x0), y1: sy(fit.slope * x0 + fit.intercept), x2: sx(x1), y2: sy(fit.slope * x1 + fit.intercept) }));
    }

    const inFamily = (m) => state.family !== null && m.families.includes(state.family);
    const order = [...mols].sort((a, b) => (inFamily(a) - inFamily(b)) || (a.bad - b.bad));
    const dots = el("svg:g", { class: state.family === null ? "dots" : "dots dimmed" });
    plotted = [];
    for (const m of order) {
      const cx = sx(axis.get(m)), cy = sy(m.logS);
      const member = inFamily(m);
      let cls = "dot";
      if (m.bad) cls += " dot--bad";
      if (member) cls += " dot--in";
      dots.append(el("svg:circle", { class: cls, cx: cx.toFixed(1), cy: cy.toFixed(1), r: member ? 4 : 3.2 }));
      coords.set(m.id, { cx, cy });
      if (state.family === null || member) plotted.push({ m, cx, cy });
    }
    svg.append(dots);
    svg.append(el("svg:g", { class: "selection", id: "selection" }));
    drawSelection();

    if (state.family !== null) {
      const f = data.families.find((x) => x.id === state.family);
      $("scatter-caption").textContent = `Highlighting ${f.n} ${f.name.toLowerCase()} molecules. Their error is ${fmt(f.rmse)} vs ${fmt(overallRmse)} overall.`;
    } else if (state.x === "logP") {
      $("scatter-caption").textContent = `Line: least-squares trend. Each +1 in logP ≈ ${sig(10 ** -fit.slope)}× less soluble (r = ${fmt(fit.r)}).`;
    } else if (state.x === "pred") {
      $("scatter-caption").textContent = "Diagonal = perfect prediction. Points above it are more soluble than the model predicted.";
    } else {
      $("scatter-caption").textContent = `Correlation with measured solubility: r = ${fmt(fit.r)}.`;
    }
  }

  function drawSelection() {
    const g = $("selection");
    const target = coords.get(state.selected);
    if (!g || !target) return;
    g.replaceChildren();
    g.append(el("svg:circle", { cx: target.cx, cy: target.cy, r: 9 }));
    const right = target.cx < svgSize().w - 140;
    g.append(el("svg:text", { x: target.cx + (right ? 13 : -13), y: target.cy - 10, "text-anchor": right ? "start" : "end" }, mols[state.selected].name));
  }

  function pointerPos(e) {
    const r = svg.getBoundingClientRect();
    return { x: e.clientX - r.left, y: e.clientY - r.top };
  }
  function nearest(pos, radius) {
    let best = null, bestD = radius * radius;
    for (const p of plotted) {
      const d = (p.cx - pos.x) ** 2 + (p.cy - pos.y) ** 2;
      if (d < bestD) { bestD = d; best = p; }
    }
    return best;
  }
  function showTooltip(p) {
    tooltip.replaceChildren(
      el("strong", {}, p.m.name),
      document.createTextNode(`measured ${fmt(p.m.logS)} · predicted ${fmt(p.m.pred)}`)
    );
    tooltip.hidden = false;
    const x = Math.min(Math.max(p.cx, 110), svgSize().w - 110);
    tooltip.style.left = `${x}px`;
    tooltip.style.top = `${p.cy}px`;
  }
  svg.addEventListener("pointermove", (e) => {
    if (e.pointerType === "touch") return;
    const p = nearest(pointerPos(e), 18);
    if (p) showTooltip(p); else tooltip.hidden = true;
  });
  svg.addEventListener("pointerleave", () => { tooltip.hidden = true; });
  svg.addEventListener("click", (e) => {
    const p = nearest(pointerPos(e), 26);
    if (p) { tooltip.hidden = true; select(p.m.id); }
  });

  /* ---------- molecule panel ---------- */
  function perLitre(logS, mw) {
    const grams = 10 ** logS * mw;
    if (grams >= 1000) return `≈ ${sig(grams)} g per litre (extremely soluble)`;
    if (grams >= 1) return `≈ ${sig(grams)} g per litre`;
    if (grams >= 1e-3) return `≈ ${sig(grams * 1e3)} mg per litre`;
    if (grams >= 1e-6) return `≈ ${sig(grams * 1e6)} µg per litre`;
    return `≈ ${sig(grams * 1e9)} ng per litre`;
  }

  function select(id) {
    state.selected = id;
    const m = mols[id];
    try { history.replaceState(null, "", `#m=${id}`); } catch { /* some file:// viewers forbid it */ }

    const structure = $("structure");
    structure.innerHTML = structures[id]; // drawn by our own build script, never user input
    structure.firstElementChild?.setAttribute("aria-label", `Structure of ${m.name}`);
    structure.firstElementChild?.setAttribute("role", "img");

    $("mol-families").replaceChildren(...m.families.map((f) => el("span", {}, data.familyNames[f])));
    $("mol-name").textContent = m.name;
    $("mol-smiles").textContent = m.smiles;

    const verdict = $("verdict");
    const box = (label, value, note) => {
      const d = el("div");
      d.append(el("span", {}, label), el("strong", {}, fmt(value)), el("small", {}, note));
      return d;
    };
    const factor = 10 ** Math.abs(m.residual);
    let sentence, cls;
    if (Math.abs(m.residual) < 0.3) { sentence = "Spot on: within 2× of the measured value."; cls = "is-good"; }
    else {
      const dir = m.residual > 0 ? "more" : "less";
      sentence = `Off by ${sig(factor)}×: it is ${dir} soluble than the model predicted.`;
      cls = m.bad ? "is-bad" : "";
    }
    verdict.replaceChildren(
      box("Measured logS", m.logS, perLitre(m.logS, m.MW)),
      box("Predicted logS", m.pred, perLitre(m.pred, m.MW)),
      el("p", { class: cls }, sentence)
    );

    const scale = $("scale");
    const lo = -12, hi = 2;
    const pct = (v) => `${((Math.min(Math.max(v, lo), hi) - lo) / (hi - lo)) * 100}%`;
    scale.replaceChildren();
    for (let t = lo; t <= hi; t += 2) {
      const tick = el("span", { class: "scale__tick" }, minus(t));
      tick.style.left = pct(t);
      scale.append(tick);
    }
    const span = el("span", { class: "scale__span" });
    span.style.left = pct(Math.min(m.logS, m.pred));
    span.style.width = `calc(${pct(Math.max(m.logS, m.pred))} - ${pct(Math.min(m.logS, m.pred))})`;
    if (!m.bad) span.style.background = "var(--muted)";
    const measured = el("span", { class: "scale__mark scale__mark--measured", title: "measured" });
    measured.style.left = pct(m.logS);
    const predicted = el("span", { class: "scale__mark scale__mark--pred", title: "predicted" });
    predicted.style.left = pct(m.pred);
    scale.append(span, predicted, measured);

    const props = [
      ["logP", fmt(m.logP)], ["Weight", `${Math.round(m.MW)} g/mol`], ["Polar surface", `${Math.round(m.TPSA)} Å²`],
      ["H-bond donors", m.HBD], ["H-bond acceptors", m.HBA], ["Aromatic rings", m.AromRings],
    ];
    $("props").replaceChildren(...props.map(([k, v]) => {
      const pair = el("div");
      pair.append(el("dt", {}, k), el("dd", {}, String(v)));
      return pair;
    }));

    const nn = mols[m.nn_index];
    const nb = $("neighbour");
    const link = el("button", { type: "button" }, nn.name);
    link.addEventListener("click", () => {
      if (state.family !== null && !nn.families.includes(state.family)) setFamily(null);
      select(nn.id);
    });
    nb.replaceChildren(document.createTextNode("Most similar molecule in the data: "), link,
      document.createTextNode(` (similarity ${m.nn_similarity.toFixed(2)})`));
    if (m.nn_similarity < 0.4) nb.append(document.createTextNode(". Nothing here is very similar, so treat this prediction with extra caution."));

    drawSelection();
  }

  /* ---------- insight bars ---------- */
  function barRow({ label, value, width, left = 0, cls = "", text, button, family }) {
    const row = el(button ? "button" : "div", { class: "bar-row" });
    if (button) { row.type = "button"; row.dataset.family = String(family); row.setAttribute("aria-pressed", "false"); }
    const track = el("span", { class: "bar-row__track" });
    const fill = el("span", { class: `bar-row__fill ${cls}`.trim() });
    fill.style.left = `${left}%`;
    fill.style.width = `${width}%`;
    const val = el("span", { class: "bar-row__value" }, text);
    if (value < 0) { val.style.right = `${100 - left}%`; } else { val.style.left = `${left + width}%`; }
    track.append(fill, val);
    row.append(el("span", { class: "bar-row__label" }, label), track);
    return row;
  }

  const lo = -1, hi = 0.45, zero = ((0 - lo) / (hi - lo)) * 100;
  const corrBars = $("corr-bars");
  corrBars.style.setProperty("--zero", `${zero}%`);
  for (const c of data.correlations) {
    const pos = ((c.r - lo) / (hi - lo)) * 100;
    corrBars.append(barRow({
      label: c.label, value: c.r, left: Math.min(pos, zero), width: Math.abs(pos - zero),
      cls: Math.abs(c.r) > 0.6 ? "is-key" : "", text: (c.r > 0 ? "+" : "") + fmt(c.r),
    }));
  }
  const ff = $("family-bars");
  const maxRmse = Math.max(...data.families.map((f) => f.rmse)) * 1.35;
  for (const f of data.families) {
    const row = barRow({
      label: f.name, value: f.rmse, width: (f.rmse / maxRmse) * 100, cls: f.rmse > overallRmse * 1.3 ? "is-hot" : "",
      text: `${fmt(f.rmse)} · n=${f.n}`, button: true, family: f.id,
    });
    row.addEventListener("click", () => {
      setFamily(state.family === f.id ? null : f.id);
      $("explorer-title").scrollIntoView({ behavior: "smooth", block: "start" });
    });
    ff.append(row);
  }
  ff.append(el("p", { class: "bars__note" }, `Average for all molecules: ${fmt(overallRmse)}`));
  const sb = $("similarity-bars");
  const maxSim = Math.max(...data.similarity.map((s) => s.rmse)) * 1.35;
  for (const s of data.similarity) {
    sb.append(barRow({
      label: `similarity ${s.bucket}`, value: s.rmse, width: (s.rmse / maxSim) * 100,
      cls: s.bucket.startsWith("<") ? "is-hot" : "", text: `${fmt(s.rmse)} · n=${s.n}`,
    }));
  }
  sb.append(el("p", { class: "bars__note" }, "Similarity = Tanimoto on Morgan fingerprints (1 = identical)."));

  /* ---------- boot ---------- */
  const fromHash = Number((location.hash.match(/m=(\d+)/) || [])[1]);
  const start = Number.isInteger(fromHash) && mols[fromHash] ? fromHash : (mols.find((m) => m.name === "Caffeine") || mols[0]).id;
  drawScatter();
  select(start);
  new ResizeObserver(() => drawScatter()).observe(frame);
})();
