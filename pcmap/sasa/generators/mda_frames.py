import ccmap

def sasa_frames_iter(md_universe, max_frame, chunk_size, step, selector, vdw_map, probe_radius=1.4, hres=False):

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
   
    if len(atom_selection.names) == 0:
        raise ValueError(f"Selector \"{selector}\" returned an empty atom set!")


    trajectory = md_universe.trajectory
    if max_frame is None:
        max_frame = len(trajectory)
    max_frame = max_frame if max_frame < len(trajectory) and max_frame > 0 else len(trajectory)
    #print(f"Processing a total of {max_frame} trajectory elements in {chunk_size} long chunks")
    steps_log = " " if step == 1 else f" [skip step {step}]"
    log = "" if not hres else "HighResolution:: "
    log += f"Computing SASA w/ a {probe_radius}A radius probe over a total of {max_frame} snapshots of {len(names)} particles each" + steps_log
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
                    resids, segids, vdw_map, probe_radius, hres )

                positions_buffer = []
        if positions_buffer:
            yield( positions_buffer, names, resnames,\
                    resids, segids, vdw_map, probe_radius, hres )                   
          
    return ( log, _iter(), max_frame )
# a thread-based function recevieving input tuple iterator and calling ccmap
def sasa_frame_task(*_args):
    args=_args[0]    
    data = ccmap.sasa_multi_mda(*args[:6], probe=args[6], hres=args[7])
    
    return data
