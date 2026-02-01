#!/usr/bin/env python
"""
Download and prepare 5-HT (serotonin) receptor PDB structures.

Serotonin targets (Purple):
- 7E2X: 5-HT1A receptor (chain R)
- 6WHA: 5-HT2A receptor (chain A)
- 4IB4: 5-HT2B receptor (chain A)
"""

from __future__ import annotations

import os
import urllib.request
from Bio.PDB import PDBParser, PDBIO, Select

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
RAW_DIR = os.path.join(BASE_DIR, "pdb_raw")
CLEAN_DIR = os.path.join(BASE_DIR, "pdb_clean")

# 5-HT receptor definitions
RECEPTORS = {
    "serotonin_5HT1A": {
        "pdb_id": "7E2X",
        "chain": "R",
        "label": "5-HT1A Receptor (7E2X)",
        "ligand": "serotonin",
        "color": "#8b5cf6",
    },
    "serotonin_5HT2A": {
        "pdb_id": "6WHA",
        "chain": "A",
        "label": "5-HT2A Receptor (6WHA)",
        "ligand": "serotonin",
        "color": "#8b5cf6",
    },
    "serotonin_5HT2B": {
        "pdb_id": "4IB4",
        "chain": "A",
        "label": "5-HT2B Receptor (4IB4)",
        "ligand": "serotonin",
        "color": "#8b5cf6",
    },
}


class ChainSelect(Select):
    """Select only specified chain and exclude water/heteroatoms."""

    def __init__(self, chain_id: str):
        self.chain_id = chain_id

    def accept_chain(self, chain):
        return chain.id == self.chain_id

    def accept_residue(self, residue):
        hetfield = residue.id[0]
        return hetfield == " "  # Standard residue only


def download_pdb(pdb_id: str, output_path: str) -> bool:
    """Download PDB file from RCSB."""
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    try:
        urllib.request.urlretrieve(url, output_path)
        return True
    except Exception as e:
        print(f"  Error downloading {pdb_id}: {e}")
        return False


def clean_pdb(input_path: str, output_path: str, chain_id: str) -> bool:
    """Extract single chain and remove water/heteroatoms."""
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure("protein", input_path)
    except Exception as e:
        print(f"  Error parsing PDB: {e}")
        return False

    io = PDBIO()
    io.set_structure(structure)
    io.save(output_path, ChainSelect(chain_id))
    return True


def main():
    os.makedirs(RAW_DIR, exist_ok=True)
    os.makedirs(CLEAN_DIR, exist_ok=True)

    print("=" * 60)
    print("Downloading and preparing 5-HT receptor structures")
    print("=" * 60)

    for name, info in RECEPTORS.items():
        pdb_id = info["pdb_id"]
        chain = info["chain"]

        print(f"\n{name}")
        print("-" * 40)
        print(f"  PDB: {pdb_id}, Chain: {chain}")
        print(f"  Label: {info['label']}")

        # Download
        raw_path = os.path.join(RAW_DIR, f"{pdb_id}.pdb")
        if not os.path.exists(raw_path):
            print(f"  Downloading {pdb_id}...")
            if not download_pdb(pdb_id, raw_path):
                continue
        else:
            print(f"  Already downloaded: {pdb_id}.pdb")

        # Clean
        clean_path = os.path.join(CLEAN_DIR, f"{name}_clean.pdb")
        print(f"  Cleaning (chain {chain})...")
        if clean_pdb(raw_path, clean_path, chain):
            size_kb = os.path.getsize(clean_path) / 1024
            print(f"  Saved: {clean_path} ({size_kb:.1f} KB)")
        else:
            print(f"  FAILED to clean {name}")

    print("\n" + "=" * 60)
    print("Done!")
    print("=" * 60)


if __name__ == "__main__":
    main()
