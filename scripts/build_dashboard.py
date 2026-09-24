"""Build the data files behind the interactive dashboard in docs/.

    python scripts/build_dashboard.py

Writes docs/data.json (numbers for every chart) and docs/structures/NNNN.svg
(one structure drawing per molecule, loaded on demand when a molecule is selected).
"""
import json
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit.Chem.Draw import rdMolDraw2D

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import chem  # noqa: E402

DOCS = ROOT / "docs"
STRUCTURES = DOCS / "structures"

DESCRIPTOR_LABELS = {
    "logP": "logP (oiliness)", "MW": "Molecular weight", "HeavyAtoms": "Heavy atoms", "AromRings": "Aromatic rings",
    "Halogens": "Halogen atoms", "AromProp": "Aromatic proportion", "RotB": "Rotatable bonds",
    "HBA": "H-bond acceptors", "TPSA": "Polar surface area", "HBD": "H-bond donors", "FracSP3": "Fraction sp3 carbon",
}
MOLECULE_COLUMNS = ["name", "smiles", "logS", "pred", "logP", "MW", "TPSA", "HBD", "HBA", "AromRings", "Halogens",
                    "nn_similarity", "nn_index", "families"]


def _short_colour(colour):
    colour = colour.upper()
    if re.fullmatch(r"#([0-9A-F])\1([0-9A-F])\2([0-9A-F])\3", colour):
        return "#" + colour[1] + colour[3] + colour[5]
    return colour


def _style(el):
    return dict(part.split(":", 1) for part in el.get("style", "").strip(";").split(";") if ":" in part)


def _minify(svg):
    """RDKit repeats a full inline style on every bond segment. Merge same-coloured strokes into
    one path and put shared defaults on a group, which shrinks each file ~5x with identical output."""
    root = ET.fromstring(svg)
    strokes, others = {}, []
    for el in root:
        tag = el.tag.split("}")[-1]
        style = _style(el)
        if tag == "rect":
            continue
        if tag == "path" and style.get("fill", "none") == "none":
            d = re.sub(r"\s*([MLQCZ])\s*", r"\1", el.get("d", "")).replace(",", " ").strip()
            width = style.get("stroke-width", "1.5px").replace("px", "")
            strokes.setdefault((_short_colour(style.get("stroke", "#000000")), width), []).append(d)
        elif tag == "path":
            others.append(f"<path d='{el.get('d')}' fill='{_short_colour(style['fill'])}' stroke='none'/>")
        elif tag == "text":
            size = style.get("font-size", "15px").replace("px", "")
            anchor = style.get("text-anchor", "start")
            extra = (f" font-size='{size}'" if size != "15" else "") + (f" text-anchor='{anchor}'" if anchor != "start" else "")
            others.append(f"<text x='{el.get('x')}' y='{el.get('y')}' fill='{_short_colour(style.get('fill', '#000000'))}'"
                          f" stroke='none'{extra}>{el.text or ''}</text>")
        else:
            attrs = " ".join(f"{k.split('}')[-1]}='{v}'" for k, v in el.attrib.items() if k != "class")
            others.append(f"<{tag} {attrs}/>")
    paths = [f"<path d='{''.join(ds)}' stroke='{colour}'" + (f" stroke-width='{w}'" if w != "1.5" else "") + "/>"
             for (colour, w), ds in strokes.items()]
    return ("<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 320 220'>"
            "<g fill='none' stroke-width='1.5' stroke-linecap='round' stroke-linejoin='round' font-family='sans-serif' font-size='15'>"
            + "".join(paths + others) + "</g></svg>")


def draw_svg(mol):
    drawer = rdMolDraw2D.MolDraw2DSVG(320, 220, -1, -1, True)  # True = plain <text> labels (much smaller files)
    opts = drawer.drawOptions()
    opts.clearBackground = False
    opts.padding = 0.08
    opts.bondLineWidth = 1.5
    opts.fixedBondLength = 26
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
    drawer.FinishDrawing()
    svg = re.sub(r"(\d+\.\d)\d+", r"\1", drawer.GetDrawingText())  # 0.1 px precision is plenty
    return _minify(svg)


def main():
    df, descriptor_cols, family_cols = chem.analyse()
    residual = df["residual"]

    def scores(pred):
        r = df["logS"] - pred
        return {"rmse": round(chem.rmse(r), 3), "r2": round(float(1 - (r ** 2).sum() / ((df.logS - df.logS.mean()) ** 2).sum()), 3)}

    families = []
    for i, name in enumerate(family_cols):
        r = residual[df[name]]
        families.append({"id": i, "name": name, "n": int(len(r)), "rmse": round(chem.rmse(r), 3),
                         "meanResidual": round(float(r.mean()), 3)})
    families.sort(key=lambda f: -f["rmse"])

    buckets = pd.cut(df["nn_similarity"], [0, 0.4, 0.6, 0.8, 1.0], labels=["< 0.4", "0.4–0.6", "0.6–0.8", "0.8–1.0"])
    similarity = [{"bucket": str(b), "n": int(len(g)), "rmse": round(chem.rmse(g), 3)}
                  for b, g in residual.groupby(buckets, observed=True)]

    corr = df[descriptor_cols + ["logS"]].corr()["logS"].drop("logS").sort_values()
    correlations = [{"descriptor": k, "label": DESCRIPTOR_LABELS[k], "r": round(float(v), 3)} for k, v in corr.items()]

    fam_matrix = df[family_cols].to_numpy()
    molecules = []
    for i, row in df.iterrows():
        molecules.append([
            row["name"], row["smiles"], round(row["logS"], 3), round(float(row["pred"]), 3), round(row["logP"], 2),
            round(row["MW"], 1), round(row["TPSA"], 1), int(row["HBD"]), int(row["HBA"]), int(row["AromRings"]),
            int(row["Halogens"]), round(float(row["nn_similarity"]), 3), int(row["nn_index"]),
            [j for j, flag in enumerate(fam_matrix[i]) if flag],
        ])

    data = {
        "source": "ESOL (Delaney, 2004) via MoleculeNet",
        "metrics": {
            "n": len(df),
            "randomForest": scores(df["pred"]),
            "esolEquation": scores(df["esol_pred"]),
            "baseline": scores(np.full(len(df), df["logS"].mean())),
            "within1": round(float((residual.abs() <= 1).mean()), 3),
        },
        "familyNames": family_cols,
        "families": families,
        "similarity": similarity,
        "correlations": correlations,
        "columns": MOLECULE_COLUMNS,
        "molecules": molecules,
    }
    DOCS.mkdir(exist_ok=True)
    (DOCS / "data.json").write_text(json.dumps(data, ensure_ascii=False, separators=(",", ":")))

    STRUCTURES.mkdir(exist_ok=True)
    for i, mol in enumerate(df["mol"]):
        (STRUCTURES / f"{i:04d}.svg").write_text(draw_svg(mol))
    print(f"Wrote docs/data.json and {len(df)} structures to docs/structures/")


if __name__ == "__main__":
    main()
