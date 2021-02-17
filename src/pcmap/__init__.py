import pypstruct, json
import ccmap as core
from .threads import run as computeMany
def contactMapThroughTransform(proteinA, proteinB,\
                               eulers, translations,\
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
    :param deserialize: ask for contact map as dictionary (True) otherwise as string.
    :type deserialize: boolean

    """
def contactMap(proteinA, proteinB=None, dist=4.5, encode=False, nThread=8, deserialize=False):
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
    try:
        dist = float(dist)
        assert dist > 0.0
    except:
        raise ValueError(f"improper distance parameter {dist}")

    try:
        if not proteinB is None:
            assert ( isinstance(proteinA, list) and isinstance(proteinB, list) )\
                or\
                   ( isinstance(proteinA, str) and isinstance(proteinB, str) )
    except:
        raise TypeError("First to pameters must have same type string or list")


    threadParam = { 'dist'   : dist, 'encode'     : False,\
                    'nThread': 8   , 'deserialize': False }

    lcMap = isinstance(proteinA, list)
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
                                      d=dist, encode=encode)
        else: # Many dimers
            threadParam.update({"pdbAtomListREC" : pdbREC,\
                                "pdbAtomListLIG" : pdbLIG})
            computeMany(**threadParam)
    else: #Monomer(s)
        if not lcMap: # Single monomer
            ccmap_as_json = core.cmap(pdbREC.atomDictorize,\
                                      d=dist, encode=encode)
        else: # Many monomers
            threadParam.update({"pdbAtomList" : pdbREC})
            computeMany(**threadParam)
    
    return json.loads(ccmap_as_json)