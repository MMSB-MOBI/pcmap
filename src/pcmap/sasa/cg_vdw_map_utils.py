import re
import numpy as np

def generate_VDW_dict_from_universe(ucg):
    vdw_d = {}
    sel = ucg.select_atoms('all')

    for res in np.unique(sel.resnames):
        R = sel.select_atoms('resname {}'.format(res))
        vdw_d[res] = get_list_names(R)
    return vdw_d

def get_dict_types(sel):
    '''Given an atom selection from MDAnalysis, 
    it will output a dictionary contaning all the beadtypes
    with their corresponding vdw radii.'''
    D = {}
    for idx, i in enumerate(sel.types):
        name = sel.names[idx]
        if re.match('^T', i):
            radii = 1.91
        elif re.match('^S', i):
            radii = 2.30
        else:
            radii = 2.64
        D[i] = radii
    return D

def get_dict_names(sel):
    '''Given an atom selection from MDAnalysis, 
    it will output a dictionary contaning all the unique beadnames
    with their corresponding vdw radii.
    NB! This is not great, if the same beadname is used with different
    beadtypes, hence different vdw radii'''
    D = {}
    for idx, i in enumerate(sel.types):
        name = sel.names[idx]
        if re.match('^T', i):
            radii = 1.91
        elif re.match('^S', i):
            radii = 2.30
        else:
            radii = 2.64
        D[name] = radii
    return D

def get_list_names(sel):
    _ = get_dict_names(sel)
    name_radii_2tup_list = [ (k,v) for k,v in _.items()]
    return name_radii_2tup_list