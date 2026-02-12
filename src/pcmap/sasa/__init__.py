import ccmap as core
from MDAnalysis import Universe

from ..utils import assert_valid_path_list, parse_selector, pdb_file_dictorize_noH
from .atom_map import atom_default_radii, validate_map
from .atom_map import atom_default_radii as vdw_AA
from .atom_map import martini3_dbg_radii as vdw_M3
from .cg_vdw_map_utils import generate_VDW_dict_from_universe
from .dto import SASA_Results
from .run_many import run as run_many

"""
Compute the Solvant accessible surface of many structures
"""


def compute_many(
    inputs,
    vdw_map=None,
    martini3=False,
    selector="all",
    npos=None,
    step=1,
    chunk_sz=5,
    probe=1.4,
    ncpu=8,
    hres=False,
):
    is_a_frame = False
    input_kwargs = {"mda_input": None, "pdb_filename_list": inputs}

    if hasattr(inputs, "universe"):
        is_a_frame = True
        if vdw_map is None:
            print("Generating VDW radii map from trajectory")
            vdw_map = generate_VDW_dict_from_universe(inputs)
        input_kwargs = {"mda_input": inputs, "pdb_filename_list": None}
    else:
        assert_valid_path_list(inputs)
        if vdw_map is None:
            if martini3:
                print("Using martini3 VDW radii map")
                vdw_map = vdw_M3
            else:
                print("Default all atoms VDW radii map")
                vdw_map = vdw_AA
    _ = run_many(
        vdw_map,
        **input_kwargs,
        max_frame=npos,
        chunk_size=chunk_sz,
        probe_radius=probe,
        ncpu=ncpu,
        selector=selector,
        step=step,
        hres=hres,
    )

    res = SASA_Results(is_a_frame)
    res.parse(_)
    return res


"""
Compute the Solvant accessible surface of one protein
"""


def compute_from_pdb(
    pdf_file_path, selector="all", martini3=False, vdw_map=None, probe=1.4, hres=False
):
    if vdw_map is None:
        if martini3:
            print("Using martini3 VDW radii map")
            vdw_map = vdw_M3
        else:
            print("Default all atoms VDW radii map")
            vdw_map = vdw_AA
    opt_segIDs = parse_selector(selector)

    noH_dict = pdb_file_dictorize_noH(pdf_file_path, opt_segIDs)
    vdw_map = validate_map(vdw_map) if not vdw_map is None else atom_default_radii
    sasa_dict = core.sasa(noH_dict, vdw_map, probe=probe, hres=hres)

    return sasa_dict
