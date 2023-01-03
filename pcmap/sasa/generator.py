import ccmap
from concurrent.futures import ThreadPoolExecutor, as_completed

from .cg_vdw_map_utils import generate_VDW_dict_from_universe
import numpy as np
from ..io import is_notebook
from tqdm import *
from tqdm.notebook import tqdm as ipython_tqdm

def sasa_frames_iter(md_universe, max_frame, chunk_size, step, selector, vdw_map, probe_radius=1.4):

    atom_selection = None
    try:      
        atom_selection = md_universe.select_atoms(selector)
    except Exception as e:
        print(f"Incorrrect atom selector \"{selector}\"")
        print(e)
        exit(1)
    
    names     = atom_selection.names
    resnames  =  atom_selection.resnames
    resids    = atom_selection.resids
    segids    = atom_selection.segids
   

    trajectory = md_universe.trajectory
    if max_frame is None:
        max_frame = len(trajectory)
    max_frame = max_frame if max_frame < len(trajectory) and max_frame > 0 else len(trajectory)
    #print(f"Processing a total of {max_frame} trajectory elements in {chunk_size} long chunks")
    steps_log = " " if step == 1 else f" [skip step { step}]"
    log = f"Processing a total of {max_frame} snapshots of {len(names)} particules each" + steps_log
    def _iter():
        positions_buffer = []
        _n = 0 # all items
        n  = 0 # items along step walk only
        for ts in trajectory[:max_frame]:
            if not _n%step == 0:
                _n += 1
                continue

            positions_buffer.append(atom_selection.positions)
            n += 1   
            if n%chunk_size == 0:
                yield( positions_buffer, names, resnames,\
                    resids, segids, vdw_map, probe_radius )

                positions_buffer = []
        if positions_buffer:
            yield( positions_buffer, names, resnames,\
                    resids, segids, vdw_map, probe_radius )                   
            #data = ccmap.np_read_multicoor(positions_buffer, names, resnames,\
            #                       resids, segids, rtype=vdw_map, probe=1.4)   
    return ( log, _iter(), max_frame )
# a thread-based function recevieving input tuple iterator and calling ccmap
def sasa_frame_task(*_args):
    args=_args[0]
    data = ccmap.np_read_multicoor(*args[:5], rtype=args[5], probe=args[6])
    
    return data

#https://stackoverflow.com/questions/51601756/use-tqdm-with-concurrent-futures
# seems ok, just need t check results order and vstack like above
def run_with_pbar(mda_universe, max_frame, chunk_size, probe_radius=1.4, ncpu=8, step=1, **kwargs):
    selector = kwargs["selector"] if "selector" in kwargs else "all"
    vdw_map = kwargs["vdw_map"] if "vdw_map" in kwargs else None

    bar_constructor = tqdm if not is_notebook() else ipython_tqdm
    # Above line is not compatible with ipywidet >= 8.0 for now ...
    bar_constructor = tqdm
    #print(f"CLI sasa Parameters: {mda_universe} {max_frame} {chunk_size} {probe_radius} {ncpu} {step}")
        
    if vdw_map is None:
        print("Generating vdw_map from trajectory...")
    else:
        raise TypeError("Tell GL he must implement a validator of VDW map")
    vdw_map = generate_VDW_dict_from_universe(mda_universe)



    res_accumul = {
        'resname' : None,
        'resID'   : None,
        'chainID' : None,
        'sasa': None
    }
    #is_notebook
    # Call generator here to get log b4 pb constructor

    (log, iter_frame, frame_num) = sasa_frames_iter(mda_universe, max_frame, chunk_size, step,\
                                    selector, vdw_map,\
                                    probe_radius=probe_radius)
    print(log)   
    with bar_constructor(total=frame_num) as pbar:
        # let's give it some more threads:
        with ThreadPoolExecutor(max_workers=ncpu) as executor:
            futures = {executor.submit(sasa_frame_task, arg): arg\
                for arg in iter_frame\
                }
            results = []
            for future in as_completed(futures):
                #arg = futures[future]
                # future seems to reduce to hash -> a number ?
                results.append( future.result() )
                pbar.update(chunk_size)
            pbar.close()
    """
        for result in results:
            for k in ['resname', 'resID', 'chainID']:
                if res_accumul[k] is None:
                    res_accumul[k] = np.array(result[k])
            for curr_pose_sasa_list in result['sasa']:
                res_accumul['sasa'] = np.array(curr_pose_sasa_list) if res_accumul['sasa'] is None\
                                      else np.vstack( [res_accumul['sasa'], curr_pose_sasa_list] )
    """
   # print(results)
    return results
    #return res_accumul
