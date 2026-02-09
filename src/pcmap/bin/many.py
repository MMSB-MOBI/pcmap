import ccmap as core
import rich_click as click
from rich import print

from ..cli_validator import (
    check_rich_and_atomic,
    path_writable,
    positive_float,
    positive_int,
)
from ..io import writeToFile
from ..threads import run as computeMany


@click.command(
    help="Compute amino acid pairwise sevelra contact over PDB structure files"
)
@click.argument("pdb_list", type=str)
@click.argument("output_file", type=str, callback=path_writable)
@click.option(
    "--dist",
    required=False,
    default=4.5,
    help=" distance threshold value for pairwise atomic contact in Anstroms",
    show_default=True,
    callback=positive_float,
)
@click.option(
    "--atomic",
    is_flag=True,
    help="output pairwise atomic contact instead of residue's",
    callback=check_rich_and_atomic,
    is_eager=True,
)
@click.option(
    "--ncpu", help="Number of cores", type=int, default=8, callback=positive_int
)
def compute_many(pdb_list: str, output_file: str, dist: float, atomic: bool, ncpu: int):
    d = {
        "nThread": ncpu,
        "deserialize": True,
        "d": dist,
        "encode": False,
        "atomic": atomic,
        "structList": pdb_list,
    }

    results = computeMany(**d)
    writeToFile(
        results,
        output_file,
    )

    print(f"[green]:heavy_check_mark: wrote {len(results)} ccmap to [i]{output_file}")
