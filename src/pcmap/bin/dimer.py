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
    help="Compute amino acid pairwise contact between two PDB structure files"
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
def compute_dimer(
    pdb_input_one: Structure,
    pdb_input_two: Structure,
    dist: float,
    rich: bool,
    atomic: bool,
):
    ccmap_as_json = core.cmap(
        pdb_input_one.atomDictorize,
        pdb_input_two.atomDictorize,
        d=dist,
        # encode=bEncode,
        atomic=atomic,
    )
    print(ccmap_as_json)


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
