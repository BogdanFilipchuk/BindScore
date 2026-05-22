#!/usr/bin/env python3
"""
estimate_protein_radius.py

Estimates the radius of a protein from a PDB file by treating it as a sphere.

Approach:
  1. Find the centroid (geometric centre) of all atoms.
  2. Compute the distance from every atom to the centroid.
  3. Return the mean of those distances as the estimated radius.

Usage:
    python estimate_protein_radius.py <file.pdb> [--chain A] [--no-hetatm] [--model 0]

Requirements:
    pip install numpy biopython
"""

import argparse
import sys
import warnings
import numpy as np

try:
    from Bio import BiopythonWarning
    from Bio.PDB import PDBParser
except ImportError:
    sys.exit(
        "BioPython is required. Install it with:\n"
        "  pip install biopython\n"
    )


def parse_args():
    p = argparse.ArgumentParser(
        description="Estimate protein radius from a PDB file (centroid + mean distance)."
    )
    p.add_argument("pdb_file", help="Path to the input PDB file")
    p.add_argument(
        "--chain", default=None,
        help="Restrict analysis to a single chain (e.g. A). Default: all chains."
    )
    p.add_argument(
        "--no-hetatm", action="store_true",
        help="Exclude HETATM records (ligands, water, etc.). Default: include them."
    )
    p.add_argument(
        "--model", type=int, default=0,
        help="Model index for NMR/multi-model files (0-based). Default: 0."
    )
    return p.parse_args()


def load_atoms(pdb_file, chain_id=None, exclude_hetatm=False, model_index=0):
    """Parse the PDB file and return an (N, 3) array of atom coordinates."""
    warnings.filterwarnings("ignore", category=BiopythonWarning)
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)

    models = list(structure.get_models())
    if model_index >= len(models):
        sys.exit(
            f"Model index {model_index} out of range "
            f"(file has {len(models)} model(s))."
        )
    model = models[model_index]

    coords = []
    for chain in model:
        if chain_id and chain.id != chain_id:
            continue
        for residue in chain:
            hetfield = residue.id[0]          # ' ' = standard, 'W' = water, 'H_*' = ligand
            if exclude_hetatm and hetfield != " ":
                continue
            for atom in residue:
                coords.append(atom.get_vector().get_array())

    if not coords:
        sys.exit(
            "No atoms found with the current filters. "
            "Check --chain and --no-hetatm options."
        )

    return np.array(coords, dtype=float)


def find_centroid(coords):
    """Return the geometric centre (mean x, y, z) of all atom positions."""
    return coords.mean(axis=0)


def estimate_radius(coords, centroid):
    """
    Compute the distance from every atom to the centroid, then return
    the mean distance as the estimated sphere radius.
    """
    distances = np.linalg.norm(coords - centroid, axis=1)
    return distances.mean(), distances


def main():
    args = parse_args()

    coords = load_atoms(
        args.pdb_file,
        chain_id=args.chain,
        exclude_hetatm=args.no_hetatm,
        model_index=args.model,
    )

    centroid = find_centroid(coords)
    radius, distances = estimate_radius(coords, centroid)

    sep = "─" * 55
    print(sep)
    print("  Protein radius estimation")
    print(sep)
    print(f"  File        : {args.pdb_file}")
    print(f"  Chain(s)    : {'all' if args.chain is None else args.chain}")
    print(f"  Atoms       : {len(coords):,}")
    print(f"  Centroid    : ({centroid[0]:.2f}, {centroid[1]:.2f}, {centroid[2]:.2f}) Å")
    print(sep)
    print(f"  Mean distance to centroid  : {radius:.2f} Å  ← estimated radius")
    print(f"  Std dev of distances       : {distances.std():.2f} Å")
    print(f"  Min distance to centroid   : {distances.min():.2f} Å")
    print(f"  Max distance to centroid   : {distances.max():.2f} Å")
    print(sep)


if __name__ == "__main__":
    main()