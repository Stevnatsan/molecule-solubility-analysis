import sys
from pathlib import Path

from rdkit import Chem

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import chem  # noqa: E402


def test_dataset_loads_every_molecule():
    df = chem.load_dataset()
    assert len(df) == 1128
    assert df["mol"].notna().all()


def test_descriptors_for_benzene():
    X = chem.featurize([Chem.MolFromSmiles("c1ccccc1")])
    assert X.loc[0, "AromRings"] == 1
    assert X.loc[0, "HeavyAtoms"] == 6
    assert X.loc[0, "AromProp"] == 1.0
    assert abs(X.loc[0, "MW"] - 78.11) < 0.05


def test_family_tags():
    mols = [Chem.MolFromSmiles(s) for s in ("CCOP(=O)(OCC)OCC", "CCO", "c1ccc2cc3ccccc3cc2c1")]
    fam = chem.families(mols, chem.featurize(mols))
    assert fam.loc[0, "Organophosphate"] and not fam.loc[1, "Organophosphate"]
    assert fam.loc[1, "Alcohol"]
    assert fam.loc[2, "Fused aromatic (3+ rings)"]


def test_nearest_neighbour_is_never_the_molecule_itself():
    mols = chem.load_dataset()["mol"].head(60)
    sim, idx = chem.nearest_neighbours(mols)
    assert all(i != j for i, j in enumerate(idx))
    assert ((sim >= 0) & (sim <= 1)).all()
