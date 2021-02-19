def isVector3F(v):
    """ Float vector of size 3 type checker """
    if not isinstance(v, list):
        return False
    if len(v) != 3:
        return False
    try:
        for x in v:
            float(x)
    except:
        return False
    return True


def isEvenEulerTranslationVectorList(e, t):
    if not ( isinstance(e, list) and isinstance(t, list) ):
        return False
    if len(e) != len(t) or not e:
        return False
    for u,v in zip(e,t):
        if not (isVector3F(u) and isVector3F(v)):
            return False
    return True