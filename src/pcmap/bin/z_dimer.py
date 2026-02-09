import json

import ccmap as core
import rich_click as click
from rich import print

from ..cli_validator import (
    Structure,
    check_rich_and_atomic,
    parse_PDB_file,
    parse_triplet,
    positive_float,
)
from ..plugins import enrich_map


@click.command(
    help="APPLY Z-DOCK like transformation to two PDB structure file and THEN compute their amino acid pairwise contacts"
)
@click.argument(
    "pdb_input_one",
    callback=parse_PDB_file,
)
@click.argument(
    "pdb_input_two",
    callback=parse_PDB_file,
)
@click.option(
    "--euler",
    help="Transformation to apply to Second PDB file, as Euler angles",
    required=True,
    callback=parse_triplet,
)
@click.option(
    "--trl",
    help="Translation to apply to Second PDB file, in cartesian space",
    required=True,
    callback=parse_triplet,
)
@click.option(
    "--ctr1",
    help="Translation to apply to CENTER First PDB structure before euler/trl transformation",
    callback=parse_triplet,
    required=True,
)
@click.option(
    "--ctr2",
    help="Translation to apply to CENTER Second PDB structure before euler/trl transformation",
    callback=parse_triplet,
    required=True,
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
def compute_z_dimer(
    pdb_input_one: Structure,
    pdb_input_two: Structure,
    euler: list[float, float, float],
    trl: list[float, float, float],
    ctr1: list[float, float, float],
    ctr2: list[float, float, float],
    dist: float,
    rich: bool,
    atomic: bool,
):
    ccmap_as_json = core.zmap(
        pdb_input_one.atomDictorize,
        pdb_input_two.atomDictorize,
        tuple(euler),
        tuple(trl),
        tuple(ctr1),
        tuple(ctr2),
        d=dist,
    )
    print(ccmap_as_json)


"""
@click.argument(
    "e2",
    help="Transformation to apply to Second PDB file, as Euler angles",
    type=float,
)
@click.argument(
    "e3",
    help="Transformation to apply to Second PDB file, as Euler angles",
    type=float,
)
@click.argument(
    "t1",
    help="Translation vector to apply to Second PDB file, in cartesian space",
    type=float,
)
@click.argument(
    "t2",
    help="Translation vector to apply to Second PDB file, in cartesian space",
    type=float,
)
@click.argument(
    "t3",
    help="Translation vector to apply to Second PDB file, in cartesian space",
    type=float,
)
@click.option(
    "--center1",
    help="Initial Translation vector to apply to first PDB file to center its coordinates ",
    callback=parse_triplet,
    default=None,
)
@click.option(
    "--center2",
    help="Initial Translation vector to apply to second PDB file to center its coordinates ",
    callback=parse_triplet,
    default=None,
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
def compute_z_dimer(
    pdb_input_one: Structure,
    pdb_input_two: Structure,
    e1,
    e2,
    e3,
    t1,
    t2,
    t3,
    center1: list[float, float, float],
    center2: list[float, float, float],
    dist: float,
    rich: bool,
    atomic: bool,
):
"""

#     ccmap_as_json = core.zmap(pdbA.atomDictorize, pdbB.atomDictorize,\
#         *vecT, **vecO, d=dist, apply=arguments['--apply'], encode=bEncode,\
#         atomic=bAtomic)


"""
if arguments['dimer']:
    if not arguments['--euler']:
        ccmap_as_json = core.cmap(pdbA.atomDictorize,\
                                  pdbB.atomDictorize,\
                                  d=dist, encode=bEncode, atomic=bAtomic)
        print(ccmap_as_json)
        exit(0)

    vecT = parseTransformVectors(arguments)
    if not '--offA' in arguments:
        ccmap_as_json = core.zmap(pdbA.atomDictorize,\
                                  pdbB.atomDictorize,\
                                  *vecT, apply=arguments['--apply'],
                                  d=dist,
                                  encode=bEncode, atomic=bAtomic)
    else:
        vecO = parseOffsetVectors(arguments)
        ccmap_as_json = core.zmap(pdbA.atomDictorize, pdbB.atomDictorize,\
            *vecT, **vecO, d=dist, apply=arguments['--apply'], encode=bEncode,\
            atomic=bAtomic)

    print(ccmap_as_json)

    if arguments['--apply']:
        # Update PDB containers
        pdbA.setCoordinateFromDictorize(pdbA.atomDictorize)
        pdbB.setCoordinateFromDictorize(pdbB.atomDictorize)
        # Dump to coordinate files
        with open("new_receptor.pdb", "w") as fp:
            fp.write( str(pdbA) )
        with open("new_ligand.pdb", "w") as fp:
            fp.write( str(pdbB) )


"""
