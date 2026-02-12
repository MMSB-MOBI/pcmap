import re
from pathlib import Path

from pypstruct import parseFilePDB

"""
    Extract segiID from a string
    Expecting "segid A or segid B ..." format
"""


def parse_selector(sel_expr):
    """
    Validate that the entire string is a sequence of segid X separated by " or ".
    Returns list all individual X's.
    """
    pattern = r"^(segid \S+)(?: or segid \S+)*$"

    if re.match(pattern, sel_expr):
        # Extract all segid individually
        matches = re.findall(r"segid \S+", sel_expr)
        return [m.split()[-1] for m in matches]

    raise ValueError(
        f'atom selector can be "all" or of the form "segid A or segid B ..." here "{sel_expr}"'
    )


"""
    Parse a pdb file and return dictorized atoms of specified segID
    without the hydrogen atoms
"""


def pdb_file_dictorize_noH(pdb_file_path, segID=None):
    pdbREC = parseFilePDB(filename=pdb_file_path)

    pdbDictREC = (
        pdbREC.chain(segID).atomDictorize if not segID is None else pdbREC.atomDictorize
    )
    noH_dict = {
        "x": [],
        "y": [],
        "z": [],
        "seqRes": [],
        "chainID": [],
        "resName": [],
        "name": [],
    }

    for x, y, z, seqRes, chainID, resName, name in zip(
        pdbDictREC["x"],
        pdbDictREC["y"],
        pdbDictREC["z"],
        pdbDictREC["seqRes"],
        pdbDictREC["chainID"],
        pdbDictREC["resName"],
        pdbDictREC["name"],
    ):
        #    if name.startswith("H"):
        #        continue
        if re.match(r"^[0-9]*H", name):
            continue
        noH_dict["x"].append(x)
        noH_dict["y"].append(y)
        noH_dict["z"].append(z)
        noH_dict["seqRes"].append(seqRes)
        noH_dict["chainID"].append(chainID)
        noH_dict["resName"].append(resName)
        noH_dict["name"].append(name)
    return noH_dict


def assert_valid_path_list(path_list):
    if not type(path_list) is list:
        raise TypeError("input for many pdb sasa is not a list")
    for _ in path_list:
        if not Path(_).is_file():
            raise TypeError(f"{_} is not a valid file")
    return True
