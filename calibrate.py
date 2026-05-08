"""
calibrate.py — Calibrate BindScore weights against experimental ΔG data
========================================================================
Reads a CSV with columns: pdb, chain_a, chain_b, dG_exp
Runs BindScore on each complex to extract interaction counts.
Fits weights by least-squares (numpy.linalg.lstsq via SVD).
Writes calibrated weights to a new interactions.json.

Usage:
    python calibrate.py --dataset pdbbind_subset.csv --output weights_calibrated.json

CSV format (no header row, or with header):
    pdb,   chain_a, chain_b, dG_exp
    1BRS,  A,       D,       -19.07
    2PTC,  E,       I,       -17.80
    ...
"""

import argparse
import csv
import json
import os
import sys
import numpy as np

# Allow running from any directory
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from BindScore import analyze, DB_PATH


FEATURE_NAMES = [
    "n_hbonds",
    "n_salt_bridges",
    "n_hydrophobic",
    "n_pi_stacks",
    "n_cation_pi",
    "buried_np_A2",
    "buried_p_A2",
]


def extract_features(result):
    return np.array([
        len(result["hbonds"]),
        len(result["salt_bridges"]),
        len(result["hydrophobic"]),
        len(result["pi_stacks"]),
        len(result["cation_pi"]),
        result["buried_np_A2"],
        result["buried_p_A2"],
    ], dtype=float)


def load_dataset(csv_path):
    rows = []
    with open(csv_path) as f:
        reader = csv.reader(f)
        for line in reader:
            if not line or line[0].startswith("#"):
                continue
            if line[0].strip().lower() == "pdb":
                continue  # header
            pdb, ca, cb, dg = [x.strip() for x in line[:4]]
            rows.append((pdb, ca, cb, float(dg)))
    return rows


def run_calibration(dataset_path, output_path, db_path=DB_PATH):
    print(f"\n  Loading dataset: {dataset_path}")
    dataset = load_dataset(dataset_path)
    print(f"  Found {len(dataset)} complexes\n")

    A_rows = []
    b_vec  = []
    failed = []

    for pdb, ca, cb, dg_exp in dataset:
        print(f"  Processing {pdb} chains {ca}/{cb} ...", end=" ", flush=True)
        try:
            result = analyze(pdb, ca, cb, db_path=db_path)
            feats  = extract_features(result)
            A_rows.append(feats)
            b_vec.append(dg_exp)
            print(f"OK  (ΔG_pred={result['dG_total']:.1f}, ΔG_exp={dg_exp:.1f})")
        except Exception as e:
            print(f"FAILED: {e}")
            failed.append(pdb)

    if len(A_rows) < 3:
        print("\n  ERROR: Need at least 3 successful complexes to calibrate.")
        sys.exit(1)

    A = np.array(A_rows)
    b = np.array(b_vec)

    # Least-squares: minimise ||Aw - b||²
    weights, residuals, rank, sv = np.linalg.lstsq(A, b, rcond=None)

    print(f"\n  Calibrated weights (kcal/mol per unit):")
    print(f"  {'Feature':<25}  {'Weight':>10}")
    print(f"  {'-'*25}  {'-'*10}")
    for name, w in zip(FEATURE_NAMES, weights):
        print(f"  {name:<25}  {w:>10.4f}")

    # Evaluate fit
    b_pred = A @ weights
    rmse   = float(np.sqrt(np.mean((b_pred - b) ** 2)))
    r2     = float(1 - np.sum((b - b_pred)**2) / np.sum((b - np.mean(b))**2))
    print(f"\n  Training RMSE : {rmse:.2f} kcal/mol")
    print(f"  R²            : {r2:.3f}")

    # Write calibrated JSON
    with open(db_path) as f:
        db = json.load(f)

    # Map weights back to interaction types
    # weights order: hbond, salt_bridge, hydrophobic, pi_stack, cation_pi,
    #                np_sasa_coeff, p_sasa_coeff
    db["hydrogen_bond"]["energy_optimal"]  = float(weights[0]) if weights[0] < 0 else db["hydrogen_bond"]["energy_optimal"]
    db["salt_bridge"]["energy_optimal"]    = float(weights[1]) if weights[1] < 0 else db["salt_bridge"]["energy_optimal"]
    db["hydrophobic"]["energy_optimal"]    = float(weights[2]) if weights[2] < 0 else db["hydrophobic"]["energy_optimal"]
    db["pi_stack"]["energy_optimal"]       = float(weights[3]) if weights[3] < 0 else db["pi_stack"]["energy_optimal"]
    db["cation_pi"]["energy_optimal"]      = float(weights[4]) if weights[4] < 0 else db["cation_pi"]["energy_optimal"]
    db["solvation"]["nonpolar_coefficient"]= float(weights[5])
    db["solvation"]["polar_coefficient"]   = float(weights[6])

    db["_calibration_meta"] = {
        "training_rmse_kcal_mol": round(rmse, 3),
        "training_r2":            round(r2, 4),
        "n_complexes":            len(A_rows),
        "failed":                 failed,
    }

    with open(output_path, "w") as f:
        json.dump(db, f, indent=2)

    print(f"\n  Calibrated database written to: {output_path}")
    if failed:
        print(f"  Warning: {len(failed)} complexes failed and were excluded: {failed}")


def build_parser():
    p = argparse.ArgumentParser(
        description="Calibrate BindScore interaction weights against experimental ΔG"
    )
    p.add_argument("--dataset", required=True,
                   help="CSV file: pdb, chain_a, chain_b, dG_exp (kcal/mol)")
    p.add_argument("--output",  default="weights_calibrated.json",
                   help="Output path for calibrated interactions.json")
    p.add_argument("--db",      default=DB_PATH,
                   help="Base interactions.json to start from")
    return p


if __name__ == "__main__":
    args = build_parser().parse_args()
    run_calibration(args.dataset, args.output, db_path=args.db)
