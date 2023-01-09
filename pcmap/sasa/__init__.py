from pypstruct import parseFilePDB
import ccmap as core
from .threads import run as compute_many_sasa
from .atom_map import atom_default_radii, validate_map
from .run_many import run as run_many
from .dto import SASA_Results
from MDAnalysis import Universe
from ..utils import pdb_file_dictorize_noH, assert_valid_path_list
from .cg_vdw_map_utils import generate_VDW_dict_from_universe
from .atom_map import martini3_dbg_radii as vdw_M3, atom_default_radii as vdw_AA
from typing import List
from pathlib import Path

def compute_many(inputs,\
                 vdw_map=None,  martini3=False,\
                 selector="all",\
                 npos=None, step=1,  
                 chunk_sz=5, probe=1.4, ncpu=8, hres=False):
    is_a_frame = False
    input_kwargs = {
        "mda_input":None,
        "pdb_filename_list" :inputs
    }
    
    if hasattr(inputs, 'universe'):
        is_a_frame = True
        if vdw_map is None:
            print("Generating VDW radii map from trajectory")
            vdw_map = generate_VDW_dict_from_universe(inputs) 
        input_kwargs = {
            "mda_input":inputs,
            "pdb_filename_list" : None
        }
    else:
        assert_valid_path_list(inputs)
        if vdw_map is None:
            if martini3:
                print("Using martini3 VDW radii map")
                vdw_map = vdw_M3
            else:
                print("Default all atoms VDW radii map")
                vdw_map = vdw_AA
    _ = run_many(vdw_map, **input_kwargs,\
            max_frame = npos, chunk_size = chunk_sz,\
            probe_radius = probe, ncpu=ncpu, selector=selector, step=step, hres=hres)
    
    res = SASA_Results(is_a_frame)
    res.parse(_)
    return res

# DEPRECATED
def compute_multi_sasa_from_frame(frame:Universe, npos=None, step=1, selector="all", vdw_map=None,\
                                    chunk_sz=5, probe=1.4, ncpu=8, hres=False):
    if vdw_map is None:
        vdw_map = generate_VDW_dict_from_universe(frame)
    else:
        vdw_map = validate_map(vdw_map)
    _ = run_many(vdw_map, mda_input=frame, max_frame = npos, chunk_size = chunk_sz,\
         probe_radius = probe, ncpu=ncpu, selector=selector, step=step, hres=hres)
    res = SASA_Results()
    res.parse(_)
    return res


"""
Compute the Solvant accessible surface of one protein
"""
def compute_from_pdb(pdf_file_path, chainID="A", vdw_map=None,probe=1.4, hres=False):
    noH_dict = pdb_file_dictorize_noH(pdf_file_path, chainID)
    vdw_map = validate_map(vdw_map) if not vdw_map is None else atom_default_radii
    sasa_dict = core.sasa(noH_dict, vdw_map,\
    probe=probe, hres=hres)
    
    return sasa_dict