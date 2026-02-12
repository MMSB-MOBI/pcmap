import pypstruct, json
import ccmap as core
from .threads import run as computeMany
from .ctype import isVector3F, isEvenEulerTranslationVectorList
from .plugins import enrich_map, assert_enrich_valid

def setThreadParameters(**kwargs):
    """ Set parameters required for parrallel computing """
    _ = setParameters(**kwargs)
    try:
        nThread = int(kwargs['nThread']) if 'nThread' in kwargs else 8
        assert nThread > 0
    except:
        raise ValueError(f"improper nThread parameter {kwargs['nThread']}")
    _['nThread'] = nThread
    _['deserialize'] = kwargs['deserialize'] if 'deserialize'\
        in kwargs else True

    return _

def setParameters(**kwargs):
    """ Set parameters required for contact map computation """
    try:
        dist = float(kwargs['dist']) if 'dist' in kwargs else 4.5
        assert dist > 0.0
    except:
        raise ValueError(f"improper distance parameter {kwargs['dist']}")
        
    param = { 
        'd'   : dist,\
        'encode' : kwargs['encode'] if 'encode' in kwargs else False,\
        'atomic' : kwargs['atomic'] if 'atomic' in kwargs else False,\
        }

    return param

def generateThroughTransform(proteinA, proteinB,\
                            euler, translation,\
                            offsetRec,\
                            offsetLig,\
                            **kwargs):
    """ Generate the PDB conformation from specified transformation"""

    if not (isVector3F(euler) and isVector3F(translation)):
        raise ValueError(f"Improper euler {euler} or translation {translation}")
    
    
    pdbRec = pypstruct.parseFilePDB(proteinA)
    pdbLig = pypstruct.parseFilePDB(proteinB)

    ccmap_as_json = core.zmap(pdbRec.atomDictorize,\
                        pdbLig.atomDictorize,\
                        euler, translation,\
                        offsetRec=offsetRec, offsetLig=offsetLig,\
                        apply=True)
    pdbRec.setCoordinateFromDictorize(pdbRec.atomDictorize)
    pdbLig.setCoordinateFromDictorize(pdbLig.atomDictorize)

    return {"receptor_oPDB": pdbRec, "ligand_oPDB": pdbLig}    

def contactMapThroughTransform(proteinA, proteinB,\
                               eulers, translations,\
                               offsetRec,\
                               offsetLig,\
                               **kwargs):
    """compute several contact map accross two provided proteins
       through the applications of provided transormations.
       Transformations are applied to the SECOND structure.

        Named parameters controls additional behaviors:
        "eulers" and "translations" allow to pass even-sized lists 
        of vectors3 which will be applied to the second structure 
        to generate dimeric conformations.
        "offsetRec" allows to pass a single translation vector to 
        center the first structure barycenter onto the origin.
        "offsetLig" allows to pass a single translation vector to 
        center the second structure barycenter onto the origin.

    :param proteinA: filepath to first structure aka "receptor"
    :type proteinA: string
    :param proteinB: filepath to second structure aka "ligand"
    :type proteinB: string
    kwargs parameters
    :param eulers: the euler angles rotation to apply to proteinB
    :type eulers: list(vector3)
    :param translations: the translation vectors to apply to proteinB
    :type eulers: list(vector3)
    
    optional
    :param nThread: thread number, default=8
    :type nThread: int
    :param deserialize: ask for contact map as dictionary (True, default) otherwise as string.
    :type deserialize: boolean

    """

    if not isEvenEulerTranslationVectorList(eulers, translations):
        raise TypeError("Transformation vector lists are of uneven length"+\
                        " and/or contain irregular values")
                        
    if not(isVector3F(offsetRec) and isVector3F(offsetLig)):
        raise TypeError("Improper offset vectors")
        

    threadParam = setThreadParameters(**kwargs)

    pdbRec = pypstruct.parseFilePDB(proteinA)
    pdbLig = pypstruct.parseFilePDB(proteinB)
        
    ccmap_as_json = computeMany(pdbA=pdbRec, pdbB=pdbLig,\
                        eulers=eulers, translations=translations,\
                        offsetRec=offsetRec, offsetLig=offsetLig,\
                        **threadParam)
    if threadParam['deserialize']:
        print("blip")
        ccmap_as_json = [ _['data'] for _ in ccmap_as_json]

    return ccmap_as_json

def contactMap(proteinA, proteinB=None, enrich=False, **kwargs):
    """compute contact map of provided structures
    
    First parameter can be PDB file OR a list of PDB files
                    
    Second parameter can be PDB file OR a list of PDB files
    
    Provided with a PDB file as single parameters:
        Compute the internal amino acid contactmap of the structure
    Provided with another PDB file as optional second parameters:
        Compute the amino acid contactmap between the two structures
    
    Provided with a list of PDB files as single parameter:
        Compute the each internal amino acid contactmap of the structure
    Provided with a list of PDB files as second parameter:
        Compute the amino acid contactmap accross pair of structures 
        at identical positions in the two lists
    
    :param proteinA : filepath to PDB files or list of path to PDB files
    :type proteinA: string|string[]
    optional
    :param proteinB : filepath to PDB files or list of path to PDB files
    :type proteinB: string|string[]
    
    :param dist: contact threshold, default=4.5 Angstroms
    :type dist: float
    :param encode: use integer encoding of contacts, default=False
    :type encode: boolean
    :param nThread: thread number, default=8
    :type nThread: int
    :param deserialize: ask for contact map as dictionary (True) otherwise as string.
    :type deserialize: boolean

    :return: contact map
    :rtype: dict
    """
    threadParam = setThreadParameters(**kwargs)
    ccmapParam  = setParameters(**kwargs)
    try:
        if not proteinB is None:
            assert ( isinstance(proteinA, list) and isinstance(proteinB, list) )\
                or\
                   ( isinstance(proteinA, str) and isinstance(proteinB, str) )
    except:
        raise TypeError("First to pameters must have same type string or list")

    lcMap = isinstance(proteinA, list)

    assert_enrich_valid(proteinA, proteinB, enrich, **kwargs)

    pdbREC = pypstruct.parseFilePDB(proteinA)\
             if not lcMap\
             else [ pypstruct.parseFilePDB(_).atomDictorize\
                    for _ in proteinA ]

    if not proteinB is None: #Dimer(s)
        pdbLIG = pypstruct.parseFilePDB(proteinB)\
                 if not lcMap\
                 else [ pypstruct.parseFilePDB(_).atomDictorize\
                        for _ in proteinB ]

        if not lcMap: # Single dimer
            ccmap_as_json = core.cmap(pdbREC.atomDictorize,\
                    pdbLIG.atomDictorize,\
                    **ccmapParam)

        else: # Many dimers
            threadParam.update({"pdbAtomListREC" : pdbREC,\
                                "pdbAtomListLIG" : pdbLIG})
            ccmap_as_json = computeMany(**threadParam)
    else: #Monomer(s)        
        if not lcMap: # Single monomer
            ccmap_as_json = core.cmap(pdbREC.atomDictorize,\
                                      **ccmapParam)
        else: # Many monomers
            threadParam.update({"pdbAtomList" : pdbREC})
            ccmap_as_json = computeMany(**threadParam)
    
    # straight calls to core lib must be deserialized
    if not lcMap and threadParam['deserialize']:
        ccmap_as_json = json.loads(ccmap_as_json)

    # Add atom informations, for students project purpose
    # We may move it to cmap package later
    if "atomic" in kwargs:
        if not lcMap and enrich and not proteinB and kwargs['atomic']:            
            return enrich_map(ccmap_as_json, pdbREC)

    return ccmap_as_json