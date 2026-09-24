"""Shared chemistry helpers for the solubility analysis.

Everything the notebook and the dashboard builder need lives here, so both
report exactly the same numbers.
"""
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdFingerprintGenerator, rdMolDescriptors
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import KFold, cross_val_predict

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[1]
DATA_PATH = ROOT / "data" / "delaney-processed.csv"

# Chemical families used to slice the errors. A molecule can belong to several.
FAMILY_SMARTS = {
    "Organophosphate": "[P](=[O,S])",
    "Nitro": "[N+](=O)[O-]",
    "Ester": "[#6]C(=O)O[#6]",
    "Amine": "[NX3;H2,H1;!$(NC=O)]",
    "Alcohol": "[CX4][OX2H]",
}


def load_dataset(path=DATA_PATH):
    """ESOL (Delaney, 2004) as distributed by MoleculeNet, with RDKit molecules attached."""
    df = pd.read_csv(path).rename(columns={
        "Compound ID": "name",
        "measured log solubility in mols per litre": "logS",
        "ESOL predicted log solubility in mols per litre": "esol_pred",
    })
    df["smiles"] = df["smiles"].str.strip()
    df["mol"] = df["smiles"].apply(Chem.MolFromSmiles)
    if df["mol"].isna().any():
        raise ValueError("Some SMILES could not be parsed")
    return df[["name", "smiles", "mol", "logS", "esol_pred"]]


def _aromatic_proportion(mol):
    return sum(atom.GetIsAromatic() for atom in mol.GetAtoms()) / mol.GetNumHeavyAtoms()


def _halogens(mol):
    return sum(atom.GetSymbol() in ("F", "Cl", "Br", "I") for atom in mol.GetAtoms())


DESCRIPTORS = {
    "logP": Crippen.MolLogP,
    "MW": Descriptors.MolWt,
    "TPSA": rdMolDescriptors.CalcTPSA,
    "HBD": Lipinski.NumHDonors,
    "HBA": Lipinski.NumHAcceptors,
    "RotB": Lipinski.NumRotatableBonds,
    "AromProp": _aromatic_proportion,
    "AromRings": rdMolDescriptors.CalcNumAromaticRings,
    "FracSP3": rdMolDescriptors.CalcFractionCSP3,
    "Halogens": _halogens,
    "HeavyAtoms": lambda m: m.GetNumHeavyAtoms(),
}


def featurize(mols):
    return pd.DataFrame({name: [fn(m) for m in mols] for name, fn in DESCRIPTORS.items()})


def families(mols, descriptors):
    """Boolean membership table: one column per chemical family."""
    out = {name: [m.HasSubstructMatch(Chem.MolFromSmarts(s)) for m in mols] for name, s in FAMILY_SMARTS.items()}
    out = pd.DataFrame(out)
    out["Heavily halogenated"] = descriptors["Halogens"].to_numpy() >= 4
    out["Fused aromatic (3+ rings)"] = descriptors["AromRings"].to_numpy() >= 3
    return out


def nearest_neighbours(mols):
    """For each molecule: Tanimoto similarity (Morgan r=2) to its closest *other* molecule, and that molecule's position."""
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fps = [gen.GetFingerprint(m) for m in mols]
    best_sim, best_idx = [], []
    for i, fp in enumerate(fps):
        sims = np.array(DataStructs.BulkTanimotoSimilarity(fp, fps))
        sims[i] = -1.0
        j = int(sims.argmax())
        best_sim.append(float(sims[j]))
        best_idx.append(j)
    return np.array(best_sim), np.array(best_idx)


def make_model(seed=0):
    return RandomForestRegressor(n_estimators=500, min_samples_leaf=2, random_state=seed, n_jobs=-1)


def cv_predictions(X, y, seed=0):
    """Out-of-fold predictions: every molecule is predicted by a model that never saw it."""
    folds = KFold(n_splits=5, shuffle=True, random_state=seed)
    return cross_val_predict(make_model(seed), X, y, cv=folds)


def rmse(residuals):
    residuals = np.asarray(residuals)
    return float(np.sqrt(np.mean(residuals ** 2)))


def analyse(path=DATA_PATH):
    """Load, featurize, predict. Returns one tidy frame used by every output."""
    df = load_dataset(path)
    X = featurize(df["mol"])
    fam = families(df["mol"], X)
    df = pd.concat([df.reset_index(drop=True), X, fam], axis=1)
    df["pred"] = cv_predictions(X, df["logS"])
    df["residual"] = df["logS"] - df["pred"]
    df["nn_similarity"], df["nn_index"] = nearest_neighbours(df["mol"])
    return df, list(X.columns), list(fam.columns)
