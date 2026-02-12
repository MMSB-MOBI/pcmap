import ccmap
from ...utils import pdb_file_dictorize_noH, parseFilePDB
import re

def pdb_list_iter(pdb_file_list, max_frame, chunk_size, step, selector, vdw_map, probe_radius=1.4, hres=False):

    opt_segIDs = None
    if selector == "all":
        opt_segIDs = None
    else :
        opt_segIDs = []
        try:
            _ = selector.split(" or ")
            for sel in _:
                opt_segIDs.append( re.match(r'^segid (\S+)', sel)[1])
        except Exception as e:
            raise ValueError(f"atom selector can be \"all\" or of the form \"segid A or segid B ...\" here \"{ selector }\"")

    
    if max_frame is None:
        max_frame = len(pdb_file_list)
    max_frame = max_frame if max_frame < len(pdb_file_list) and max_frame > 0 else len(pdb_file_list)
    #print(f"Processing a total of {max_frame} trajectory elements in {chunk_size} long chunks")
    steps_log = " " if step == 1 else f" [skip step {step}]"
    log = "" if not hres else "HighResolution:: "
    log += f"Computing SASA w/ a {probe_radius}A radius probe over a total of {max_frame} pdb records " + steps_log
    def _iter():
        dictorized_atoms_buffer = []
        _n = 0 # all items
        n  = 0 # items along step walk only
        for pdb_file in pdb_file_list[:max_frame]:
            if not _n%step == 0:
                _n += 1
                continue
            dictorized_atoms_buffer.append( pdb_file_dictorize_noH(pdb_file, opt_segIDs) )
            if len(dictorized_atoms_buffer['name']) == 0:
                raise ValueError(f"Selector \"{selector}\" returned an empty atom set!")

            n += 1   
            if n%chunk_size == 0:
                yield( dictorized_atoms_buffer, vdw_map, probe_radius, hres )

                dictorized_atoms_buffer = []
        if dictorized_atoms_buffer:
            yield( dictorized_atoms_buffer, vdw_map, probe_radius, hres )                   
          
    return ( log, _iter(), max_frame )
# a thread-based function recevieving input tuple iterator and calling ccmap
def sasa_pdb_list_task(*_args):
    args=_args[0]
    data = ccmap.sasa_list(args[0], args[1], probe=args[2], hres=args[3])
    
    return data
