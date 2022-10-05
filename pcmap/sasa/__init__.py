from pypstruct import parseFilePDB
import ccmap as core
from .threads import run as compute_many_sasa
from .atom_map import atom_default_radii


def setThreadParameters(**kwargs):
    """ Set parameters required for sasa computation  """
    try:
        probeRadius = float(kwargs['r']) if 'r' in kwargs else 1.4
        assert probeRadius > 0.0
    except:
        raise ValueError(f"improper probe radius parameter {kwargs['r']}")
        
    return { 
        'r'   : probeRadius,
        'atom_radii' : kwargs['rtype'] if 'rtype' in kwargs else atom_default_radii
        }

# Basic single thread monomeric sasa computation

### We need hydro filter out ?

def compute_many(protein, **kwargs):
    """compute solvant accessible surface of provided structures
    
    First parameter can be PDB file OR a list of PDB files
                    
    
    Provided with a PDB file as single parameters:
        Compute the solvant accessible surface of the structure
    Provided with a list of PDB files as single parameter:
        Compute the solvant accessible surface of each structure
    
    :param protein: filepath to PDB files or list of path to PDB files
    :type protein: string|string[]
    
    :param nThread: thread number, default=8
    :type nThread: int
    
    :return: a list of dictionary
    :rtype: list
    """
   # threadParam = setThreadParameters(**kwargs)
    sasaBaseThreadParam  = setThreadParameters(**kwargs)
    try:
        assert ( isinstance(protein, list) or isinstance(protein, str) )
    except:
        raise TypeError("The \"protein\" pameter must have same type string or list")

    input = None
    #many structure
    if ( isinstance(protein, list) ):
        
        try :
            input = [ parseFilePDB(_).atomDictorize\
                    for _ in protein ]
        except:
            raise TypeError("Error while parsing pdb files")
        results = compute_many_sasa(\
            input, **sasaBaseThreadParam, **kwargs)
        return results
    # single structure
    try :
        input =  parseFilePDB(protein).atomDictorize
    except:
        raise TypeError("Error while parsing the pdb file")
    
    sasa_dict = core.sasa(input,\
        sasaBaseThreadParam['atom_radii'],
        probe=sasaBaseThreadParam['r']
        ) # convert to array
    return sasa_dict  

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