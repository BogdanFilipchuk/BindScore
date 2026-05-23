"""
protein_solvation_entropy.py
────────────────────────────────────────────────────────────────────────────────
Calculates the entropy and thermodynamic solvation quantities for a protein
(solute) in water (solvent), following the framework of:

    Sun, Q. (2022). "The Hydrophobic Effects: Our Current Understanding."
    Molecules, 27(20), 7009.  https://doi.org/10.3390/molecules27207009
    PMC9609269

Physics summary
───────────────
When a protein is dissolved in water an interface forms between the two phases.
The paper partitions water into:
  • Interfacial water  – topmost layer at the protein surface; tends to form
                         weaker DA (single-donor/single-acceptor) H-bonds.
  • Bulk water         – retains tetrahedral DDAA H-bonds (stronger, lower S).

Hydration free energy (Eq. 6):
    ΔG_hyd = ΔG_water-water + 8·ΔG_DDAA·r_H₂O / R_protein

A critical radius Rc (Eq. 7) divides two solvation regimes:
  • R < Rc  → initial solvation  : interfacial term dominates → solutes disperse
  • R > Rc  → hydrophobic solvation: bulk term dominates      → solutes aggregate

Water constants are taken directly from the paper (293 K, 0.1 MPa):
  ΔH_DDAA = 11.35 kJ/mol   ΔS_DDAA = 29.66 J/(mol·K)   (van't Hoff fit, Fig. 5)
  ΔG_water-water = −1500 cal/mol = −6276 J/mol           (ref [85] in paper)

The protein is modelled as a sphere whose radius is derived from its molecular
mass (computed from the amino-acid sequence in the PDB file) and the average
protein density (1350 kg/m³).

Usage
─────
    python protein_solvation_entropy.py myprotein.pdb
    python protein_solvation_entropy.py myprotein.pdb --temp 310
    python protein_solvation_entropy.py myprotein.pdb --temp 310 --chain A
"""

import math
import sys
from pathlib import Path
sys.path.append("src/bindscore/pdb_file_treatment")
import bindscore.pdb_file_treatment.protein_radius_estimation as radius_estimation 
#aaaaaaaaaaaaaaaaaaaaaaaaaaaa
try:
    from Bio.PDB import PDBParser, PPBuilder
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
except ImportError:
    sys.exit(
        "ERROR: Biopython is required.\n"
        "Install it with:  pip install biopython"
    )
# ══════════════════════════════════════════════════════════════════════════════
#  Fixed water constants  (Sun 2022, all at 293 K / 0.1 MPa unless noted)
# ══════════════════════════════════════════════════════════════════════════════
R_GAS          = 8.314          # J mol⁻¹ K⁻¹
CAL_TO_J       = 4.184          # 1 cal = 4.184 J

# Van't Hoff fit to Raman OH intensities (paper Eq. 3 / Fig. 5)
# These are for the DDAA (tetrahedral) → free-OH transition in pure water.
DH_DDAA        = 11.35e3        # J/mol
DS_DDAA        = 29.66          # J/(mol·K)

# Gibbs free energy of bulk water (ref [85] of the paper)
DG_WATER_WATER = -1500 * CAL_TO_J   # J/mol  =  −6276 J/mol

# Effective radius of a single water molecule  (paper Sec 3: d = 3.8 Å → r = 1.9 Å)
R_H2O          = 1.9e-10        # m

# Average number of H-bonds per molecule in the DDAA network
N_HB           = 2

# Average protein density (well-established literature value)
PROTEIN_DENSITY = 1350.0        # kg/m³

# Avogadro constant
AVOGADRO       = 6.02214076e23


# ══════════════════════════════════════════════════════════════════════════════
#  PDB → protein radius
# ══════════════════════════════════════════════════════════════════════════════

def parse_protein(pdb_path: str, chain_id: str | None = None) -> tuple[str, str, float]:
    """
    Parse a PDB file and return (protein_name, sequence, sphere_radius_m).

    The sphere-equivalent radius R is derived from the molecular weight M via:
        V  = M / (ρ · Nₐ)          [volume of one molecule]
        R  = (3V / 4π)^(1/3)

    Parameters
    ----------
    pdb_path : path to the .pdb file
    chain_id : if given, only residues from that chain are used

    Returns
    -------
    name     : PDB structure ID
    sequence : one-letter amino-acid sequence
    radius   : sphere-equivalent radius (metres)
    """
    path = Path(pdb_path)
    if not path.exists():
        sys.exit(f"ERROR: file not found: {pdb_path}")

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(path.stem, str(path))
    name = path.stem.upper()

    # Build peptide chains
    ppb = PPBuilder()
    sequence = ""
    for model in structure:
        for pp in ppb.build_peptides(model, aa_only=True):
            # Filter by chain if requested
            if chain_id is not None:
                chain = pp[0].get_parent().id
                if chain != chain_id:
                    continue
            sequence += str(pp.get_sequence())

    if not sequence:
        sys.exit(
            f"ERROR: no amino-acid residues found in {pdb_path}"
            + (f" for chain {chain_id}" if chain_id else "")
            + ".\nCheck the file or try --chain with a different letter."
        )

    # Molecular weight (Da = g/mol) from sequence
    analysis = ProteinAnalysis(sequence)
    mw_gmol  = analysis.molecular_weight()   # g/mol
    mw_kgmol = mw_gmol * 1e-3               # kg/mol

    coords = radius_estimation.load_atoms(pdb_path,chain_id=chain_id)
    centroid = radius_estimation.find_centroid(coords)
    radius = radius_estimation.estimate_radius(coords, centroid)[0]
    return name, sequence, radius


# ══════════════════════════════════════════════════════════════════════════════
#  Water thermodynamics at a given temperature
# ══════════════════════════════════════════════════════════════════════════════

def dG_DDAA(T: float) -> float:
    """
    Gibbs free energy of the DDAA tetrahedral H-bond at temperature T (K).
    ΔG = ΔH − T·ΔS   (negative = stabilising)
    """
    return DH_DDAA - T * DS_DDAA   # J/mol


def critical_radius(T: float) -> float:
    """
    Critical solute radius Rc (m) — paper Eq. 7.
    Rc = 8·|ΔG_DDAA|·r_H₂O / |ΔG_water-water|
    """
    return 8.0 * abs(dG_DDAA(T)) * R_H2O / abs(DG_WATER_WATER)


# ══════════════════════════════════════════════════════════════════════════════
#  Thermodynamic quantities for the protein–water system
# ══════════════════════════════════════════════════════════════════════════════

def interfacial_to_volume_ratio(R: float) -> float:
    """
    Ratio of interfacial water molecules to total water volume for a sphere
    of radius R (m):   ratio = 4·r_H₂O / R   (paper Sec 3)
    """
    return 4.0 * R_H2O / (R*1e-10)


#def dG_solute_water(R: float, T: float) -> float:
    """
    Gibbs free energy of the solute-water interface (J/mol) — Eq. 5.
    ΔG_sw = ΔG_DDAA(T) · (4·r_H₂O/R) · N_HB
    """
    return dG_DDAA(T) * interfacial_to_volume_ratio(R) * N_HB


#def dG_hydration(R: float, T: float) -> float:
    """
    Total hydration free energy (J/mol) — Eq. 6.
    ΔG_hyd = ΔG_water-water + 8·ΔG_DDAA(T)·r_H₂O / R
    """
    return DG_WATER_WATER + 8.0 * dG_DDAA(T) * R_H2O / R


def dS_interfacial(R: float, T: float) -> float:
    """
    Entropy contribution of interfacial water (J mol⁻¹ K⁻¹).

    DA bonds at the interface are weaker (lower ΔH) than DDAA bulk bonds,
    so they carry higher configurational entropy.  The gain relative to bulk
    scales with the interfacial fraction:
        ΔS_int = (ΔH_DDAA / T) · ratio
    """
    return (DH_DDAA / T) * interfacial_to_volume_ratio(R)


def dS_bulk(R: float) -> float:
    """
    Entropy contribution of bulk water around the protein (J mol⁻¹ K⁻¹).

    Represents the DDAA ordering cost in the bulk fraction:
        ΔS_bulk = ΔS_DDAA · (1 − ratio)
    """
    ratio = interfacial_to_volume_ratio(R)
    return DS_DDAA * (1.0 - ratio)


#def dS_total(R: float, T: float) -> float:
    """
    Net entropy change of water upon protein solvation (J mol⁻¹ K⁻¹).
    ΔS_total = ΔS_interfacial − ΔS_bulk
    """
    return dS_interfacial(R, T) - dS_bulk(R)


#def dH_solvation(R: float, T: float) -> float:
    """
    Enthalpy change of solvation (J/mol), from ΔG = ΔH − T·ΔS:
        ΔH = ΔG_solute-water + T·ΔS_interfacial
    """
    return dG_solute_water(R, T) + T * dS_interfacial(R, T)


def solvation_regime(R: float, T: float) -> str:
    """
    Classify the solvation regime (paper Eq. 7 / Sec 3):
      initial (R < Rc)    → interfacial term dominates → protein stays dispersed
      hydrophobic (R ≥ Rc) → bulk term dominates       → tendency to aggregate
    """
    Rc = critical_radius(T)
    return "initial / dispersed (R < Rc)" if R < Rc else "hydrophobic / aggregated (R ≥ Rc)"


# ══════════════════════════════════════════════════════════════════════════════
#  Output
# ══════════════════════════════════════════════════════════════════════════════

def print_results(name: str, sequence: str, radius_m: float, T: float) -> None:
    R   = radius_m
    Rc  = critical_radius(T)
    #dG  = dG_DDAA(T)

    # Pre-compute all quantities
    ratio   = interfacial_to_volume_ratio(R)
    #dGsw    = dG_solute_water(R, T)
    #dGhyd   = dG_hydration(R, T)
    dSint   = dS_interfacial(R, T)
    dSblk   = dS_bulk(R)
    #dStot   = dS_total(R, T)
    #dHsolv  = dH_solvation(R, T)
    #neg_TdS = -T * dStot
    regime  = solvation_regime(R, T)

    w = 52   # column width

    def row(label, value, unit=""):
        print(f"  {label:<{w-2}} {value}  {unit}")

    print()
    print("═" * w)
    print(f"  Protein–Water Solvation Entropy  ·  {name}")
    print("═" * w)

    print("\n── Input ──────────────────────────────────────────")
    row("PDB / protein name",         name)
    row("Sequence length",            f"{len(sequence)} residues")
    row("Temperature",                f"{T:.1f}", "K")

    print("\n── Protein geometry ───────────────────────────────")
    row("Sphere-equivalent radius R", f"{R:.2f}", "Å")
    row("Sphere-equivalent radius R", f"{R * 0.1:.3f}",  "nm")

    print("\n── Water constants  (Sun 2022, 0.1 MPa) ───────────")
    #row("ΔH_DDAA  (van't Hoff, Raman)", f"{DH_DDAA / 1e3:.2f}", "kJ/mol")
    row("ΔS_DDAA  (van't Hoff, Raman)", f"{DS_DDAA:.2f}",       "J/(mol·K)")
    #row("ΔG_DDAA  at given T",          f"{dG / 1e3:.4f}",      "kJ/mol")
    #row("ΔG_water-water",               f"{DG_WATER_WATER/1e3:.4f}", "kJ/mol")
    row("Critical radius Rc",           f"{Rc * 1e10:.2f}", "Å")

    print("\n── Solvation geometry ─────────────────────────────")
    row("Interfacial/volume ratio (4·r_H₂O/R)", f"{ratio:.4f}")
    row("Solvation regime",            regime)

    print("\n── Entropy decomposition ──────────────────────────")
    row("ΔS_interfacial  (DA bond gain)",  f"{dSint:+.4f}", "J/(mol·K)")
    row("ΔS_bulk         (DDAA order cost)", f"{dSblk:+.4f}", "J/(mol·K)")
    #row("ΔS_total        (net, water)",    f"{dStot:+.4f}", "J/(mol·K)")
    #row("−T·ΔS_total",                    f"{neg_TdS/1e3:+.4f}", "kJ/mol")

    #print("\n── Free energy & enthalpy ─────────────────────────")
    #row("ΔG_solute-water  (Eq. 5)",   f"{dGsw  / 1e3:+.4f}", "kJ/mol")
    #row("ΔG_hydration     (Eq. 6)",   f"{dGhyd / 1e3:+.4f}", "kJ/mol")
    #row("ΔH_solvation",              f"{dHsolv/ 1e3:+.4f}", "kJ/mol")
    #row("ΔG = ΔH − T·ΔS  (check)",   f"{(dHsolv + neg_TdS)/1e3:+.4f}", "kJ/mol")

    print()
    print("═" * w)
    print()


