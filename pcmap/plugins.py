def enrich_map(ccmap_in, pdbObj):
    """ Add x,y,z coordinates to atoms of provided atomic contact map"""    
    d = {}
    av = pdbObj.atomVectorize
    trace = [ (_.seqRes, _.chainID) for _ in pdbObj.trace ]

    for x, y, z, seqRes, chainID, resName, name in zip(*av[:]):
        d[ f"{name}{resName}{seqRes}{chainID}" ] = (name, resName, seqRes, chainID, x, y, z)

    ccmap_out = {"type" : "atomic_rich", "data": []}
    for a1, a2, dist in ccmap_in["data"]:
        #print (a1, a2)
        #print(d[''.join(a1)],  d[''.join(a2)], dist)
        ccmap_out["data"]. append( (d[''.join(a1)],  d[''.join(a2)], dist) )
    
    return ccmap_out

def assert_enrich_valid(proteinA, proteinB, enrich, **kwargs):
    if not enrich:
        return True

    if proteinB:
        raise ValueError("enrich contact map cannot be computed on two bodies "\
                       + f'(proteinB is not None "{proteinB}"")')
    
    atomic_error_msg = 'Enriching contact map requires atomic contact map computation, '\
                     + 'please set "atomic" named parameter to True'
    if not "atomic" in kwargs:
        raise ValueError(atomic_error_msg)
    
    if not kwargs["atomic"]:
        raise ValueError(atomic_error_msg)
    
    return True