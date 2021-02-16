def tripletParser(value):
    try :
        v = [ float(_) for _ in value.split(',') ]
    except:
        print(f"Cannot parse vector of 3 real value from {value}")
        return None
    return v

def transformationFileParser(fileName):
    pass


def parseTransformVectors(parameters):
    vecT = []
    for vType in ['--euler', '--trans']:
        try:
            vecT.append( tripletParser(parameters[vType]) )
        except:
            raise ValueError (f"Can't parse {parameters[vType]} as the {vType} vector")
    return vecT

def parseOffsetVectors(parameters):
    vecO = {}
    for vType, vKey in zip(['--offA', '--offB'], ['offsetRec', 'offsetLig']):
        try:
            vecO[vKey] = tripletParser(parameters[vType])
        except:
            raise ValueError(f"Can't parse {parameters[vType]} as the {vType} vector")
    print(vecO)
    return vecO
    