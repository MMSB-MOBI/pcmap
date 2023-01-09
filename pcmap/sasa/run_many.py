import ccmap
from concurrent.futures import ThreadPoolExecutor, as_completed


import numpy as np
from ..io import is_notebook
from .generators.mda_frames import sasa_frames_iter, sasa_frame_task
from .generators.pdb_lists import pdb_list_iter, sasa_pdb_list_task
from tqdm import *
#from tqdm.notebook import tqdm as ipython_tqdm

#from tqdm.autonotebook import tqdm

def run(vdw_map, mda_input=None, pdb_filename_list=None,\
        max_frame=None, chunk_size=None, \
        selector = None,\
        probe_radius=1.4, ncpu=8, step=1, hres=False):
   # Defautls settings to PDB inputs
    generator = pdb_list_iter 
    task      = sasa_pdb_list_task
    iter_args = (pdb_filename_list, max_frame, chunk_size, step,\
                    selector, vdw_map)
    if (mda_input is None and pdb_filename_list is None) or (not mda_input is None and not pdb_filename_list is None):
        raise ValueError("Please provide a mda Universe object or a list of pdb file path")
    
    iter_kwargs = {"probe_radius":probe_radius, "hres":hres}
    if not mda_input is None:
        generator = sasa_frames_iter
        task      = sasa_frame_task
        iter_args = (mda_input, max_frame, chunk_size, step,\
                                    selector, vdw_map)
    else: # pdb_file_path_list
        if not vdw_map:
            raise ValueError("Please provide many pdb calculation with a vdw radii map")

    results = run_with_pbar(generator, task, iter_args, iter_kwargs, ncpu, chunk_size)
    return results

#https://stackoverflow.com/questions/51601756/use-tqdm-with-concurrent-futures
# seems ok, just need t check results order and vstack like above
def run_with_pbar(iter,  task, iter_args, iter_kwargs, n_worker, chunk_sz):


    bar_constructor = tqdm

    #is_notebook
    # Call generator here to get log b4 pb constructor

    (log, iter_frame, frame_num) = iter(*iter_args, **iter_kwargs)
    print(log)

    inc_step = chunk_sz if chunk_sz < frame_num else frame_num
    with bar_constructor(total=frame_num) as pbar:
        # let's give it some more threads:
        with ThreadPoolExecutor(max_workers=n_worker) as executor:
            futures = {executor.submit(task, arg): arg\
                for arg in iter_frame\
                }
            results = []
            for future in as_completed(futures):
                #arg = futures[future]
                # future seems to reduce to hash -> a number ?
                results.append( future.result() )
                pbar.update(inc_step)
            pbar.close()
    return results
