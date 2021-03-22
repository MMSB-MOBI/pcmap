"""Compute amino acid contact map within a single protein or across two proteins
Usage:
    pcmap single <proteinA> [--distance=<distance> --encode]
    pcmap dimer  <proteinA> <proteinB> [--distance=<distance> --encode]  
    pcmap dimer  <proteinA> <proteinB> --euler=<euler_triplet> --trans=<translation_triplet> [(--offA=<offsetA> --offB=<offsetB>)] [--distance=<distance> --encode --apply]     
    pcmap many   (--structures=<structureList> | <proteinA> <proteinB> <transformation_file>) [--distance=<distance> --ncpu=<thread_num> --output=<filename> --encode]
    pcmap -h | --help

Options:
  -h --help     Show this screen.
  proteinA: protein structure in PDB format
  proteinB: protein structure in PDB format
  euler_triplet: three comma separated values specifting the vector of Euler angle rotation to apply to proteinB
  translation_triplet: three comma separated values specifying the translation vector to apply to proteinB
  offsetA: translation vector centering proteinA barycenter, offsetB must be provided.
  offsetB: translation vector centering proteinB barycenter, offsetA must be provided.
  --distance: distance threshold value for pairwise atomic contact, default=4.5 Angstroms
  --apply: apply provided tansformation to proteinA and proteinB and write their coordinates
  --cpu: thread number, default=8
  --encode: encode amino acid contact as integers, default=False
  --output: many contact map file output, default="contact_map_many.json"
"""

from docopt import docopt
import pypstruct
from .io import * 
import ccmap as core
from .threads import run as computeMany

arguments = docopt(__doc__)
#print(arguments)

dist = 4.5
if arguments['--distance']:
    try:
        dist = float(arguments['--distance'])
        assert dist > 0.0
    except:
        print (f"Can't parse {arguments['--distance']} as contact threshold distance")
        exit(1)

pdbA = None
if not arguments['--structures']:
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

bEncode = arguments['--encode']
if arguments['many']:
    try:
        nThread = int(arguments['--ncpu']) if arguments['--ncpu'] else 8
        d = {"nThread": nThread, "deserialize": True,\
             "d":dist, "encode": bEncode}
    except:
        print("--npu arguments is not an integer")
        exit(1)
    if arguments['--structures']:
        d.update( {"structList" : arguments['--structures']} )
    else:
        d.update({"pdbA" : pdbA,
                  "pdbB" : pdbB, 
                  "transformation" : arguments['<transformation_file>'] 
                })
    results = computeMany(**d)
    writeToFile(results, arguments['--output'])
    exit(0)

if arguments['single']:
    ccmap_as_json = core.cmap(pdbA.atomDictorize, d=dist, encode=bEncode)
    print(ccmap_as_json)
    exit(0)

if arguments['dimer']:
    if not arguments['--euler']:
        ccmap_as_json = core.cmap(pdbA.atomDictorize,\
                                  pdbB.atomDictorize,\
                                  d=dist, encode=bEncode)
        print(ccmap_as_json)
        exit(0)    
    
    vecT = parseTransformVectors(arguments)
    if not '--offA' in arguments:
        ccmap_as_json = core.zmap(pdbA.atomDictorize,\
                                  pdbB.atomDictorize,\
                                  *vecT, apply=arguments['--apply'],
                                  d=dist,
                                  encode=bEncode)
    else:  
        vecO = parseOffsetVectors(arguments)
        ccmap_as_json = core.zmap(pdbA.atomDictorize, pdbB.atomDictorize,\
            *vecT, **vecO, d=dist, apply=arguments['--apply'], encode=bEncode)

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