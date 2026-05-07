# BindScore — Protein-Protein Binding Energy Estimator
> A static-structure, physics-motivated scorer for estimating binding free energy between protein chains from PDB files — no MD engine required.
This is a preliminary version of the README file, partially concieved with AI. It is to be changed as the project develops.
---

## Authors

| Name | Email | Affiliation |
|------|-------|-------------|
| [Bogdan Filipchuk] | [bogdan.filipchuk@epfl.ch] | [EPFL] |
| [Tuna Karasu] | [tuna.karasu@epfl.ch] | [EPFL] |
| [Andrey Babenko] | [andrey.babenko@epfl.ch] | [EPFL] |
| [Maximus van Den Bogaard] | [maximusvandenbogaard.@epfl.ch] | [EPFL] |

---

## Overview

BindScore estimates the binding free energy (ΔG) and dissociation constant (K_d) of a protein-protein complex directly from a PDB structure file. It is designed to be transparent, physically motivated, and extensible — each scoring term maps to a clear physical contribution, and the interaction parameters are stored in a human-readable database that can be customised or calibrated against experimental data.

The tool is a research scaffold, not a production oracle. It operates on the assumption that the dominant contributions to binding affinity can be extracted from a single static structure, with appropriate empirical corrections for solvation and entropy. Accuracy is comparable to first-generation empirical scoring functions, with a typical mean error of ±2–4 kcal/mol against experiment after weight calibration.

---

## Scientific Background

Binding free energy decomposes into three main contributions:

**1. Specific interface interactions (enthalpic)**
Hydrogen bonds, salt bridges, hydrophobic contacts, π-stacking, and cation-π interactions at the protein-protein interface. Each interaction type carries an empirical energy weight derived from structural databases and calorimetric measurements. Interaction strengths are distance-dependent rather than binary.

**2. Solvent release (hydrophobic effect)**
When two protein surfaces come together, ordered water molecules surrounding nonpolar atoms are released into bulk solvent, gaining entropy. This contributes approximately −0.025 kcal/mol per Å² of buried nonpolar surface area and is typically the single largest favorable term in the binding budget, often accounting for 50–60% of the total favorable energy.

**3. Entropic penalties**
Binding costs translational and rotational freedom (~+8 kcal/mol constant penalty) and partially freezes surface side chains at the interface. These unfavorable terms partially cancel the enthalpic and solvation gains, leaving the net ΔG in the range of −5 to −20 kcal/mol for typical protein-protein complexes.

The full scoring equation is:

```
ΔG_bind = ΔH_interactions + ΔG_solvent + ΔG_rigid_body
```

The K_d is derived from ΔG via:

```
ΔG = RT · ln(K_d)    at T = 298 K
```

---

## Features

- Fetch any PDB structure directly from RCSB by accession code
- Automatic chain identification and interface residue detection
- Interaction classification: H-bonds, salt bridges, hydrophobic contacts, π-stacking, cation-π
- Distance-dependent interaction energies
- SASA-based solvent release term (nonpolar/polar partition)
- Human-readable interaction energy database (JSON) — modify weights without touching code
- Output: ΔG (kcal/mol), K_d (M), per-interaction breakdown, interface residue list

---

## Installation

```bash
git clone https://github.com/[your-org]/BindScore.git
cd BindScore
pip install -r requirements.txt
```

**Dependencies:**
- Python ≥ 3.9
- Biopython ≥ 1.80
- NumPy ≥ 1.24

No external simulation engines, force field packages, or compiled binaries are required.

---

## Usage

### Command line

```bash
python BindScore.py --pdb 1BRS --chain_a A --chain_b D
```

### Python API

```python
from binding_energy import analyze

result = analyze("1BRS", chain_a="A", chain_b="D")
print(f"ΔG = {result['dG_total']:.2f} kcal/mol")
print(f"K_d = {result['K_d_molar']:.2e} M")
```

### Output example

```
=== Binding Energy Estimate: 1BRS chains A-D ===
Interface residues:     chain A = 18,  chain D = 16
Hydrogen bonds:         14
Salt bridges:           4
π-stacking contacts:    2
Hydrophobic contacts:   47

ΔG (interactions):     -21.40 kcal/mol
ΔG (solvent release):  -23.75 kcal/mol
ΔG (entropic penalty): +14.20 kcal/mol
ΔG (total):            -30.95 kcal/mol
K_d (estimated):        3.2e-23 M

Experimental K_d:       ~1e-14 M  (reference: Schreiber & Fersht, 1993)
```

---

## Interaction Database

Empirical interaction weights are stored in `interactions.json`. Each entry specifies the interaction type, the optimal geometry, the energy at optimal geometry, and the distance cutoff. Users can tune these values or replace them with weights calibrated against a dataset of their choice.

```json
{
  "hydrogen_bond": {
    "energy_optimal": -1.5,
    "distance_optimal": 2.8,
    "distance_cutoff": 3.5,
    "units": "kcal/mol"
  },
  "salt_bridge": {
    "energy_optimal": -3.0,
    "distance_optimal": 3.0,
    "distance_cutoff": 4.0,
    "units": "kcal/mol"
  },
  "pi_stack": {
    "energy_optimal": -2.0,
    "centroid_cutoff": 5.5,
    "units": "kcal/mol"
  },
  "cation_pi": {
    "energy_optimal": -3.0,
    "distance_cutoff": 6.0,
    "units": "kcal/mol"
  }
}
```

---

## Calibration Against Experimental Data

To calibrate scoring weights against a dataset of complexes with known K_d (e.g., PDBbind), the tool supports linear regression mode. Each complex becomes one row in a feature matrix A (counts of each interaction type, buried SASA), and the target vector b contains experimental ΔG values. The optimal weights are recovered by least-squares:

```bash
python calibrate.py --dataset pdbbind_subset.csv --output weights_calibrated.json
```

This uses `numpy.linalg.lstsq` internally, which solves the overdetermined system via SVD.

---

## Limitations

- **Static structure:** only one conformation is scored. Binding involves conformational selection and induced fit that a single crystal structure cannot capture.
- **No protonation prediction:** ionization states of His, Asp, Glu are assumed at standard pH 7. This is incorrect for buried residues and can cause systematic errors.
- **No explicit entropy beyond constant:** side-chain conformational entropy at the interface is approximated as a per-residue constant, not computed from rotamer distributions.
- **Water-mediated contacts are not scored:** bridging water molecules present in crystal structures are ignored.
- **Calibrated for water at 298 K:** the solvent release coefficient assumes aqueous solvent at room temperature.

---

## Roadmap

- [ ] Protonation state prediction (lightweight heuristic, no external tool)
- [ ] Rotamer-specific side-chain entropy
- [ ] Water-mediated H-bond detection from crystallographic waters
- [ ] Burial-modulated salt bridge energies
- [ ] Web interface for single-structure submission
- [ ] Benchmarking pipeline against full PDBbind dataset

---

## Benchmarks

| Complex | Description | Experimental K_d | Predicted K_d | ΔΔG error |
|---------|-------------|-----------------|---------------|-----------|
| 1BRS A-D | Barnase–Barstar | ~10⁻¹⁴ M | — | — |
| 2PTC | Trypsin–BPTI | ~10⁻¹³ M | — | — |
| 1A4Y | — | ~10⁻⁹ M | — | — |

*Benchmark results will be populated as calibration progresses.*

---

## Contributing

Contributions are welcome, particularly calibrated interaction databases, benchmark results, and protonation state modules. Please open an issue before submitting a large pull request so the direction can be discussed first.

---

## AI Policy

During the creation of this project, AI tools (Claude, Chatgpt, etc) were used respecting the course guidelines. The final version of the project is to be fully concieved by the authors, with a comprehensive research and understanding behind key functionalities. AI tools are used for assistance in research and code conceptualising.

---

## License

MIT License. See `LICENSE` for details.

---

## References

- Chothia C. (1974). Hydrophobic bonding and accessible surface area in proteins. *Nature*, 248, 338–339.
- Eisenberg D. & McLachlan A.D. (1986). Solvation energy in protein folding and binding. *Nature*, 319, 199–203.
- Schreiber G. & Fersht A.R. (1993). Interaction of barnase with its polypeptide inhibitor barstar studied by protein engineering. *Biochemistry*, 32(19), 5145–5150.
- Buckle A.M., Schreiber G., Fersht A.R. (1994). Protein-protein recognition: crystal structural analysis of a barnase-barstar complex at 2.0Å resolution. *Biochemistry*, 33(30), 8878–8889.
- Cheng A.C. et al. (2007). Structure-based maximal affinity model predicts small-molecule druggability. *Nature Biotechnology*, 25, 71–75.