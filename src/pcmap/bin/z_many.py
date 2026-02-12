import ccmap as core
import rich_click as click
from rich import print

from ..cli_validator import (
    Structure,
    check_rich_and_atomic,
    is_json,
    parse_PDB_file,
    path_writable,
    positive_float,
    positive_int,
)
from ..io import writeToFile
from ..threads import run as computeMany


@click.command(
    help="Perform several transformation of references structure and compute their amino acid pairwise contacts"
)
@click.argument(
    "receptor_PDB_file", type=click.Path(exists=True), callback=parse_PDB_file
)
@click.argument(
    "ligand_PDB_file", type=click.Path(exists=True), callback=parse_PDB_file
)
@click.argument("transformation_file", type=click.Path(exists=True), callback=is_json)
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
def compute_many(
    receptor_pdb_file: Structure,
    ligand_pdb_file: Structure,
    transformation_file: str,
    output_file: str,
    dist: float,
    atomic: bool,
    ncpu: int,
):
    d = {
        "nThread": ncpu,
        "deserialize": True,
        "d": dist,
        "encode": False,
        "atomic": atomic,
        "transformation": transformation_file,
        "pdbA": receptor_pdb_file,
        "pdbB": ligand_pdb_file,
    }

    results = computeMany(**d)
    writeToFile(
        results,
        output_file,
    )

    print(f"[green]:heavy_check_mark: wrote {len(results)} ccmap to [i]{output_file}")
