#!/usr/bin/env python
"""Create biotin and serotonin 3D structures as SDF files."""

import os
from rdkit import Chem
from rdkit.Chem import AllChem

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

LIGANDS = {
    "biotin": {
        "smiles": "C1[C@H]2[C@@H]([C@@H](S1)CCCCC(=O)O)NC(=O)N2",
        "color": "#ec4899",  # Pink
    },
    "serotonin": {
        "smiles": "NCCc1c[nH]c2ccc(O)cc12",
        "color": "#8b5cf6",  # Purple
    },
}


def generate_3d_sdf(smiles: str, name: str, output_path: str) -> bool:
    """Generate 3D conformer and save as SDF."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"  Error: Invalid SMILES for {name}")
        return False

    mol = Chem.AddHs(mol)

    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    result = AllChem.EmbedMolecule(mol, params)
    if result == -1:
        print(f"  Warning: Could not embed {name}, trying with random coords")
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())

    try:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
    except Exception:
        print(f"  Warning: MMFF optimization failed for {name}")

    mol.SetProp("_Name", name)

    writer = Chem.SDWriter(output_path)
    writer.write(mol)
    writer.close()

    return True


def main():
    print("Creating biotin and serotonin ligand files...")

    for name, info in LIGANDS.items():
        smiles = info["smiles"]
        output_path = os.path.join(BASE_DIR, f"{name}.sdf")

        print(f"\n{name.upper()}")
        print(f"  SMILES: {smiles}")

        if generate_3d_sdf(smiles, name, output_path):
            print(f"  Saved: {output_path}")
        else:
            print(f"  FAILED to generate {name}")

    print("\nDone!")


if __name__ == "__main__":
    main()
