"""
cli.py
======
Command-line interface for BindScore.

Entry point registered in pyproject.toml as the `bindscore` command.
"""

from __future__ import annotations
import click


# ---------------------------------------------------------------------------
# Root group
# ---------------------------------------------------------------------------

@click.group()
@click.version_option(package_name="bindscore")
def cli():
    """BindScore — fast estimation of protein-protein binding entropy."""



# ---------------------------------------------------------------------------
# `bindscore getchains` - command for getting chains 
# ---------------------------------------------------------------------------


@cli.command()
@click.argument("pdb_ID")
def getchains(pdb_id:str):
    """
    Prints the chains available in the PDB.
    """
    import pathlib
    from bindscore.pdb_file_treatment.pdb_utils_fetch import fetch_pdb_file, fetch_pdb_data
    from bindscore.parsing.pdb_utils_protein import Protein_Structure

    pdb_path = pathlib.Path(pdb_id)

    if not pdb_path.exists():
        click.echo(f"Fetching {pdb_id.upper()} from RCSB...")
        pdb_path = fetch_pdb_file(pdb_id)

    myprotein = Protein_Structure(fetch_pdb_data(str(pdb_path)))
    click.echo(f"Polypeptide chains in {myprotein.get_ID()}: {sorted(myprotein.chains())}")
    if(len(myprotein.small_molecules())!=0):
        click.echo(f"Small molecules in {myprotein.get_ID()} : {sorted(myprotein.small_molecules())}")






# ---------------------------------------------------------------------------
# `bindscore score` - main scoring command
# ---------------------------------------------------------------------------

@cli.command()
@click.argument("pdb",        metavar="PDB")
@click.argument("chain_a",    metavar="CHAIN_A")
@click.argument("chain_b",    metavar="CHAIN_B")
@click.option("--temperature", "-T", default=300.0, show_default=True,
              help="Temperature in Kelvin.")
@click.option("--breakdown",   "-b", is_flag=True, default=False,
              help="Print per-component breakdown.")
@click.option("--json",        "as_json", is_flag=True, default=False,
              help="Output results as JSON.")
@click.option("--show_failed",     "-sf",      is_flag=True, default=False,
              help="Show detailed error messages from failed submodules.")
def score(pdb, chain_a, chain_b, temperature, breakdown, as_json, show_failed):
    """Compute the interaction deltaG components for 2 chains of a PDB. You already need to know the chains so run getchains first

    PDB      Path to a PDB file, or a 4-character RCSB ID (auto-fetched).
    CHAIN_A  Chain ID of the first binding partner  (e.g. A).
    CHAIN_B  Chain ID of the second binding partner (e.g. B).

    Example:

        bindscore score 1BRS.pdb A D --breakdown
    """
    import pathlib
    from bindscore.pdb_file_treatment.pdb_utils_fetch import fetch_pdb_file
    from bindscore.scoring import compute_total_entropy

    # Resolve PDB: local file or remote fetch
    pdb_path = pathlib.Path(pdb)
    if not pdb_path.exists():
        click.echo(f"Fetching {pdb.upper()} from RCSB...")
        pdb_path = fetch_pdb_file(pdb)

    click.echo(f"Scoring {pdb_path.name}  chains {chain_a}/{chain_b}  T={temperature} K")

    result = compute_total_entropy(
        str(pdb_path), chain_a, chain_b,
        T=temperature,
        return_breakdown=True,
    )

    minus_T_dS_kcal = -(result.dS_total * temperature) / 4184.0

    if result.failed_modules:
        if show_failed:
            for msg in result.failed_modules:
                click.echo(f"  [warn] {msg}")
        else:
            n = len(result.failed_modules)
            click.echo(f"  Warning: {n} submodule(s) failed and contributed 0 (use -v for details)")

    if as_json:
        import json
        data = {
            "note":                  "preliminary - work in progress",
            "dS_total_J_per_molK":   result.dS_total,
            "dS_trans_rot":          result.dS_trans_rot,
            "dS_hydrophobic":        result.dS_hydrophobic,
            "dS_sidechain":          result.dS_sidechain,
            "dS_backbone":           result.dS_backbone,
            "minus_T_dS_kcal_mol":   minus_T_dS_kcal,
        }
        click.echo(json.dumps(data, indent=2))
        return

    # ── Delta-H section ──────────────────────────────────────────────────────
    dH_kcal = None
    try:
        import sys, importlib.util
        from pathlib import Path

        # The parsing and scoring modules still use bare-name imports
        # (e.g. `from pdb_utils_protein import *`), so we add their
        # directories to sys.path before importing them.
        pkg_dir = Path(importlib.util.find_spec("bindscore").origin).parent
        for directory in [pkg_dir / "parsing", pkg_dir / "scoring"]:
            if str(directory) not in sys.path:
                sys.path.insert(0, str(directory))

        from bindscore.pdb_file_treatment.pdb_utils_fetch import fetch_pdb_data
        from bindscore.parsing.pdb_utils_protein import Protein_Structure
        from bindscore.parsing.pdb_utils_inter import Interaction
        from bindscore.scoring.pdb_utils_enthalpy import binding_enthalpy

        pdb_data    = fetch_pdb_data(str(pdb_path))
        protein     = Protein_Structure(pdb_data)
        inter       = Interaction(protein, chain_a, chain_b, threshold=5.0)
        breakdown_h = binding_enthalpy(inter)
        dH_kcal     = breakdown_h.get("TOTAL", 0.0)

        if breakdown:
            for k, v in breakdown_h.items():
                if k == "TOTAL":
                    continue
                click.echo(f"  {k:20s}  {v:+12.4f}  kJ/mol")
            click.echo("  " + "-" * 46)


            
        click.echo("")
        click.echo("  Delta-H  (preliminary - work in progress)")
        click.echo("  " + "-" * 46)
        click.echo(f"  {'TOTAL dH':20s}  {dH_kcal:+12.4f} kJ/mol")
    except Exception as exc:
        click.echo(f"  dH not available: {exc}")



     # ── Delta-S section ──────────────────────────────────────────────────────
    click.echo("")
    click.echo("  Delta-S  (preliminary - work in progress)")
    click.echo("  " + "-" * 46)
    if breakdown:
        click.echo(f"  {'Trans+Rot':20s}  {result.dS_trans_rot:+12.2f}  J/(mol*K)")
        if result.trans_rot_detail:
            d = result.trans_rot_detail
            click.echo(f"    trans={d.translational:+.2f}  rot={d.rotational:+.2f}")

        click.echo(f"  {'Hydrophobic':20s}  {result.dS_hydrophobic:+12.2f}  J/(mol*K)")
        if result.hydrophobic_detail:
            d = result.hydrophobic_detail
            click.echo(f"    R_a={d.R_a:.1f} A  R_b={d.R_b:.1f} A")

        click.echo(f"  {'Sidechain':20s}  {result.dS_sidechain:+12.2f}  J/(mol*K)")
        if result.sidechain_detail:
            click.echo(f"    n_interface={result.sidechain_detail.n_interface_residues}")

        click.echo(f"  {'Backbone (NMA)':20s}  {result.dS_backbone:+12.2f}  J/(mol*K)")
        if result.backbone_detail:
            click.echo(f"    n_modes={result.backbone_detail.n_modes_matched}")

        click.echo("  " + "-" * 46)

    click.echo(f"  {'TOTAL dS':20s}  {result.dS_total:+12.2f}  J/(mol*K)")
    click.echo(f"  {'-T*dS':20s}  {minus_T_dS_kcal:+12.4f}  kcal/mol")




# ---------------------------------------------------------------------------
# `bindscore web` - launch the Streamlit UI !
# ---------------------------------------------------------------------------

@cli.command()
@click.option("--port", "-p", default=8501, show_default=True,   #HONESTLY have no idea, chat recommended to add it, to figure out later
              help="Port to serve the app on.")
def app(port):
    """Launch the BindScore Streamlit web interface."""
    import subprocess, importlib.util
    from pathlib import Path

    pkg_directory = Path(importlib.util.find_spec("bindscore").origin).parent
    app_path = pkg_directory / "web" / "app.py"

    if not app_path.exists():
        raise click.ClickException(f"Web app not found at {app_path}")

    click.echo(f"Starting BindScore web UI on port {port}...")
    subprocess.run(["streamlit", "run", str(app_path), "--server.port", str(port)])


# ---------------------------------------------------------------------------
# `bindscore fetch` - download a PDB file
# ---------------------------------------------------------------------------

@cli.command()
@click.argument("pdb_id", metavar="PDB_ID")
@click.option("--out", "-o", default=".", show_default=True,
              help="Directory to save the file.")
def fetch(pdb_id, out):
    """Download a PDB file from RCSB and save it locally.

    Example:

        bindscore fetch 1BRS --out ./pdbs/
    """
    import pathlib
    from bindscore.pdb_file_treatment.pdb_utils_fetch import fetch_pdb_file

    save_dir = pathlib.Path(out)
    save_dir.mkdir(parents=True, exist_ok=True)
    path = fetch_pdb_file(pdb_id, save_dir=save_dir)
    click.echo(f"Saved -> {path}")
