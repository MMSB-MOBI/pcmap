import json

import ccmap as core
import rich_click as click
from rich import print

from ..cli_validator import (
    Structure,
    check_rich_and_atomic,
    parse_PDB_file,
    positive_float,
)
from ..plugins import enrich_map


@click.command(
    help="Compute amino acid pairwise contact in a single chain PDB structure file"
)
@click.argument(
    "pdb_input",
    callback=parse_PDB_file,
)
@click.option(
    "--dist",
    required=False,
    default=4.5,
    help=" distance threshold value for pairwise atomic contact in Anstroms",
    show_default=True,
    callback=positive_float,
)
@click.option(
    "--rich",
    is_flag=True,
    help="add cartesian coordinates to contact map output",
    is_eager=True,
    callback=check_rich_and_atomic,
)
@click.option(
    "--atomic",
    is_flag=True,
    help="output pairwise atomic contact instead of residue's",
    callback=check_rich_and_atomic,
    is_eager=True,
)
def compute_monomer(pdb_input: Structure, dist: float, rich: bool, atomic: bool):
    ccmap_as_json = core.cmap(
        pdb_input.atomDictorize, d=dist, encode=False, atomic=atomic
    )
    if rich:
        _ = json.loads(ccmap_as_json)
        ccmap_as_json = json.dumps(enrich_map(_, pdb_input), indent=4)

    print(ccmap_as_json)
