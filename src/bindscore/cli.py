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
# `bindscore score` — main scoring command
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
def score(pdb, chain_a, chain_b, temperature, breakdown, as_json):
    """Compute binding entropy ΔS for a PDB complex.

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
        click.echo(f"Fetching {pdb.upper()} from RCSB…")
        pdb_path = fetch_pdb_file(pdb)

    click.echo(f"Scoring {pdb_path.name}  chains {chain_a}/{chain_b}  T={temperature} K")

    result = compute_total_entropy(
        str(pdb_path), chain_a, chain_b,
        T=temperature,
        return_breakdown=True,
    )

    if as_json:
        import json
        data = {
            "dS_total":       result.dS_total,
            "dS_trans_rot":   result.dS_trans_rot,
            "dS_hydrophobic": result.dS_hydrophobic,
            "dS_sidechain":   result.dS_sidechain,
            "dS_backbone":    result.dS_backbone,
        }
        click.echo(json.dumps(data, indent=2))
        return

    click.echo("")
    if breakdown:
        click.echo(f"  {'Trans+Rot':20s}  {result.dS_trans_rot:+12.2f}  J/(mol·K)")
        if result.trans_rot_detail:
            d = result.trans_rot_detail
            click.echo(f"    trans={d.translational:+.2f}  rot={d.rotational:+.2f}")

        click.echo(f"  {'Hydrophobic':20s}  {result.dS_hydrophobic:+12.2f}  J/(mol·K)")
        if result.hydrophobic_detail:
            d = result.hydrophobic_detail
            click.echo(f"    R_a={d.R_a:.1f} Å  R_b={d.R_b:.1f} Å")

        click.echo(f"  {'Sidechain':20s}  {result.dS_sidechain:+12.2f}  J/(mol·K)")
        if result.sidechain_detail:
            click.echo(f"    n_interface={result.sidechain_detail.n_interface_residues}")

        click.echo(f"  {'Backbone (NMA)':20s}  {result.dS_backbone:+12.2f}  J/(mol·K)")
        if result.backbone_detail:
            click.echo(f"    n_modes={result.backbone_detail.n_modes_matched}")

        click.echo(f"  {'─'*20}  {'─'*12}")

    click.echo(f"  {'TOTAL ΔS':20s}  {result.dS_total:+12.2f}  J/(mol·K)")
    minus_T_dS_kcal = -(result.dS_total * temperature) / 4184.0
    click.echo(f"  {'-T·ΔS':20s}  {minus_T_dS_kcal:+12.4f}  kcal/mol")


# ---------------------------------------------------------------------------
# `bindscore fetch` — download a PDB file
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
    click.echo(f"Saved → {path}")
