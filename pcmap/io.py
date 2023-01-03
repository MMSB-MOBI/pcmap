import json
import pypstruct

def is_notebook() -> bool:
    try:
        shell = get_ipython().__class__.__name__
        if shell == 'ZMQInteractiveShell':
            return True   # Jupyter notebook or qtconsole
        elif shell == 'TerminalInteractiveShell':
            return False  # Terminal running IPython
        else:
            return False  # Other type (?)
    except NameError:
        return False      # Probably standard Python interpreter

def writeToFile(data, fname):
    """Dump data object as JSON formated string in optional filename, 
    default value for file is 'contact_map_many.json'"""
    with open(fname if fname else 'contact_map_many.json', "w") as fp:
        json.dump({"contact maps" : data}, fp, indent = 6)
    
def tripletParser(value):
    """ Parse comma separated string of numbers in a vector3(float)"""
    try :
        v = [ float(_) for _ in value.split(',') ]
    except:
        print(f"Cannot parse vector of 3 real value from {value}")
        return None
    return v

def parseTransformVectors(parameters):
    """ Parse Euler and translation vectors
    
    :param file: parameters dict
    :type file: dict
    :return: dict(eulerTriplet, transTriplet)
    :rtype: list(vector3, vector3)
    """
    vecT = []
    for vType in ['--euler', '--trans']:
        try:
            vecT.append( tripletParser(parameters[vType]) )
        except:
            raise ValueError (f"Can't parse {parameters[vType]} as the {vType} vector")
    return vecT

def parseOffsetVectors(parameters):
    """ Parse offset vectors
    
    :param file: parameters dict
    :type file: dict
    :return: dict()
    :rtype: dict('offsetRec':vector3, 'offsetLig':vector3)
    """
    vecO = {}
    for vType, vKey in zip(['--offA', '--offB'], ['offsetRec', 'offsetLig']):
        try:
            vecO[vKey] = tripletParser(parameters[vType])
        except:
            raise ValueError(f"Can't parse {parameters[vType]} as the {vType} vector")
    return vecO

def parseTransformationFile(file):
    """ Parse a JSON coordinates transformation file
    
    :param file: path to JSON file
    :type file: string
    :return: tuple()
    :rtype: tuple(list(vector3), list(vector3), vector3, vector3)
    """
    with open(file, "r") as fp:
        vectors = json.load(fp)
        eulers     = [tuple(_) for _ in vectors['euler']]
        translations  = [tuple(_) for _ in vectors['translation']]
        ligOffset = vectors['ligOffset']
        recOffset = vectors['recOffset'] 
    return eulers, translations, recOffset, ligOffset

def parseStructFileTabList(file):
    """ Parse a tabulated list of PDB structure fileName
    
    :param file: path to tabulated file
    :type file: string
    :return: tuple()
    :rtype: tuple(int, list([atomVectors, atomVectors]|[atomVectors]))
    """
    structList = []
    fieldCount = None
    with open(file, "r") as fp:
        for l in fp:
            _ = l.rstrip().split()
            if not ( len(_) == 1 or len(_) == 2):
                raise ValueError(f"Irregular structure names at {_} on input {file}")
            if fieldCount is None:
                fieldCount = len(_)
            if fieldCount != len(_):
                raise ValueError(f"Irregular structure count (got{len(_)} expected {fieldCount}) at {_} on input {file}")
            structList.append(
                        [ pypstruct.parseFilePDB(x).atomDictorize for x in _ ]
                        )
    return fieldCount, structList

