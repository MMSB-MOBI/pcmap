from pypstruct import parseFilePDB
import ccmap as core
from .atom_map import atom_default_radii

# Basic single thread monomeric sasa computation

### To edit for sasa

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

###

def compute_from_pdb(pdf_file_path):
    pdb_container     = parseFilePDB(filename=pdf_file_path)
    pdb_as_atom_dict  = pdb_container.atomDictorize
    noH_dict = {
        "x" : [],
        "y" : [],
        "z" : [],
        "seqRes" : [],
                        "chainID" : [],
                        "resName" : [],
                        "name" : []
    }

    for x, y, z, seqRes, chainID, resName, name in zip(\
        pdb_as_atom_dict["x"], pdb_as_atom_dict["y"], pdb_as_atom_dict["z"],\
        pdb_as_atom_dict["seqRes"], pdb_as_atom_dict["chainID"],\
        pdb_as_atom_dict["resName"], pdb_as_atom_dict["name"]
        ):
        if name.startswith("H"):
            continue
        noH_dict["x"].append(x)
        noH_dict["y"].append(y)
        noH_dict["z"].append(z)
        noH_dict["seqRes"].append(seqRes)
        noH_dict["chainID"].append(chainID)
        noH_dict["resName"].append(resName)
        noH_dict["name"].append(name)

    sasa_dict = core.sasa(noH_dict, atom_default_radii)
    
    return sasa_dict