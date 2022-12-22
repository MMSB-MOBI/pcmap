"""Compute amino acid contact map within a single protein or across two proteins
Usage:
    pcmap single <proteinA> [--distance=<distance> --encode --atomic --rich]
    pcmap dimer  <proteinA> <proteinB> [--distance=<distance> --encode --atomic]  
    pcmap dimer  <proteinA> <proteinB> --euler=<euler_triplet> --trans=<translation_triplet> [(--offA=<offsetA> --offB=<offsetB>)] [--distance=<distance> --encode --apply --atomic]     
    pcmap many   (--structures=<structureList> | <proteinA> <proteinB> <transformation_file>) [--distance=<distance> --ncpu=<thread_num> --output=<filename> --encode --atomic]
    pcmap sasa <tpr_file> <xtc_file> [--npos=<number> --csz=<number>]
    pcmap -h | --help

Options:
  -h --help     Show this screen.
  proteinA: protein structure in PDB format
  proteinB: protein structure in PDB format
  euler_triplet: three commas separated values specifting the vector of Euler angle rotation to apply to proteinB
  translation_triplet: three commas separated values specifying the translation vector to apply to proteinB
  offsetA: translation vector centering proteinA barycenter, offsetB must be provided.
  offsetB: translation vector centering proteinB barycenter, offsetA must be provided.
  --distance: distance threshold value for pairwise atomic contact, default=4.5 Angstroms
  --apply: apply provided tansformation to proteinA and proteinB and write their coordinates
  --cpu: thread number, default=8
  --encode: encode amino acid contacts as integers, default=False
  --rich: add cartesian coordinates to contact map output, only compatible with one single body atomic computation, default=False
  --atomic: output pairwise atomic contact instead of residue's, default=False
  --output: many contact map or many sasa file output, default="contact_map_many.json" OR "sasa_many.json"
  TPR_FILE: simulation file TPR
  XTC_FILE: simulation file XTC
  --npos: number of snapshot to process, default=10 
  --csz : chunk size, default=5

"""

from docopt import docopt
from pcmap.plugins import assert_enrich_valid
import pypstruct
from .io import * 
from .sasa.generator import run_with_pbar as compute_sasa_frame
import ccmap as core
from .threads import run as computeMany
from .plugins import assert_enrich_valid, enrich_map
import json
import MDAnalysis as md

arguments = docopt(__doc__)
#print(arguments)

if arguments['sasa']:
    n_poses = 10
    if arguments['--npos']:
        try:
            n_poses = int(arguments['--npos'])
            assert n_poses > 0
        except:
            print (f"Can't parse {arguments['--npos']} as number of poses")
            exit(1)

    chunk_sz = 5
    if arguments['--csz']:
        try:
            chunk_sz = int(arguments['--csz'])
            assert chunk_sz > 0
        except:
            print (f"Can't parse {arguments['--csz']} as number of poses per chunck")
            exit(1)
    print("Reading trajectory files...")
    ucg = md.Universe(arguments['<tpr_file>'],arguments['<xtc_file>'])
    data = compute_sasa_frame(ucg, n_poses, chunk_sz, probe_radius=1.4, ncpu=8)
    fout = arguments['--output'] if arguments['--output'] else "sasa_many.json"
    with open(fout, 'w') as fp:
        fp.write(str(data))
    exit(0)

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
bAtomic = arguments['--atomic']
bEnrich = arguments['--rich']

if arguments['many']:
    try:
        nThread = int(arguments['--ncpu']) if arguments['--ncpu'] else 8
        d = {"nThread": nThread, "deserialize": True,\
             "d":dist, "encode": bEncode, "atomic" : bAtomic}
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
    if bEnrich and not bAtomic:
        raise ValueError('Enriching contact map requires atomic contact map computation, '\
                       + 'please set "--atomic" flag')
    ccmap_as_json = core.cmap(pdbA.atomDictorize,\
                            d=dist, encode=bEncode, atomic=bAtomic)
    if bEnrich:
        _ = json.loads(ccmap_as_json)
        ccmap_as_json = json.dumps( enrich_map(_, pdbA), indent=4)
    print(ccmap_as_json)
    exit(0)

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