"""Compute amino acid contact map within a single protein or across two proteins
Usage:
    pcmap single <proteinA> [--distance]
    pcmap dimer  <proteinA> <proteinB> [--distance --apply]  
    pcmap dimer  <proteinA> <proteinB> --euler=<euler_triplet> --trans=<translation_triplet> [(--offA=<offsetA> --offB=<offsetB>)] [--distance --apply]     
    pcmap many   (--structures=<structureList> | <proteinA> <proteinB> <transformation_file>) [--distance]

Options:
  -h --help     Show this screen.
  proteinA: protein structure in PDB format
  proteinB: protein structure in PDB format
  euler_triplet: three comma separated values specifting the vector of Euler angle rotation to apply to proteinB
  translation_triplet: three comma separated values specifying the translation vector to apply to proteinB
  offsetA: translation vector centering proteinA barycenter, offsetB must be provided.
  offsetB: translation vector centering proteinB barycenter, offsetA must be provided.
  --distance  : distance threshold value for pairwise atomic contact, default=4.5 Angstroms
  --apply  : apply provided tansformation to proteinA and proteinB and write their coordinates
"""

from docopt import docopt
import pypstruct
from .io import * 
import ccmap as core

arguments = docopt(__doc__)
print(arguments)

dist = 4.5
if arguments['--distance']:
    try:
        dist = float(arguments['--distance'])
        assert dist > 0.0
    except:
        print (f"Can't parse {arguments['--distance']} as contact threshold distance")
        exit(1)

pdbA = None
try:
    pdbA = pypstruct.parseFilePDB(arguments['<proteinA>'])
    assert not pdbA is None
except:
    print (f"Can't parse {arguments['<proteinA>']} as first PDB record")
    exit(1)

pdbB = None
if arguments['<proteinB>']:
    try:
        pdbB = pypstruct.parseFilePDB(arguments['<proteinB>'])
        assert not pdbB is None
    except:
        print (f"Can't parse {arguments['<proteinB>']} as second PDB record")
        exit(1)


if arguments['many']:
    

if arguments['single']:
    ccamp_as_json = core.cmap(pdbA.atomDictorize, d=dist, encode=False)
    print(ccamp_as_json)
    exit(0)

if arguments['dimer']:
    if not arguments['--euler']:
        ccamp_as_json = core.cmap(pdbA.atomDictorize, pdbB.atomDictorize, d=dist, encode=False)
        print(ccamp_as_json)
        exit(0)    
    
    vecT = parseTransformVectors(arguments)
    if not '--offA' in arguments:
        ccamp_as_json = core.zmap(pdbA.atomDictorize, pdbB.atomDictorize , *vecT, apply=arguments['--apply'] )
    else:  
        vecO = parseOffsetVectors(arguments)
        ccamp_as_json = core.zmap(pdbA.atomDictorize, pdbB.atomDictorize , *vecT, **vecO, apply=arguments['--apply'] )

    print(ccamp_as_json)

    if arguments['--apply']:
        # Update PDB containers from previous examples
        pdbA.setCoordinateFromDictorize(pdbA.atomDictorize)
        pdbB.setCoordinateFromDictorize(pdbB.atomDictorize)
        # Dump to coordinate files
        with open("new_receptor.pdb", "w") as fp:
            fp.write( str(pdbA) )
        with open("new_ligand.pdb", "w") as fp:
            fp.write( str(pdbB) )