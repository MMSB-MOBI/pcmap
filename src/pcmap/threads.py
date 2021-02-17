import sys, threading, json, time, copy, os, ccmap
import io
from .io import parseStructFileTabList, parseTransformationFile
from .generators import *

def isTransformAsFile(parameters):
    """ Test if desired operation relies of application 
        of transformation provided as a file
    """
    return ('transformation' in parameters)

def isTransformAsVector(parameters):
    """ Test if desired operation relies on applications 
        of transformations specified as euler, translation 
        and offset vectors
    """
    return ('eulers' in parameters and 'translation' in parameters and\
        'offsetRec' in parameters and 'offsetLig' in parameters )

def isExplicitAsFile(parameters):
    """ Test if desired operation relies on straight 
        computations over a list of PDB provided in a file
    """
    return 'fStructList' in parameters
           
def isExplicitAsVector(parameters):
    """ Test if desired operation relies on straight 
        computations over a list of PDB provided as atom vectors
    """
    return ('pdbAtomListREC' in parameters and\
            'pdbAtomListLIG' in parameters) or\
           ('pdbAtomList' in parameters)


def isExplicitAsMonomersVector(parameters):
    return 'pdbAtomList' in parameters

def isExplicitAsDimersVector(parameters):
    return 'pdbAtomList' in parameters

def run(threadNum=8, **kwargs):
    """ Run threads to process supplied structures 
    
        Returns a list of dictionaries, one per input structure/transformation
        The dictionaries hold the computed contact maps, input order is preserved.
        Inside each dictionary the value corresponding to the "data" key 
        is the entier structure contact map. It stores a list of small dictionaries 
        holding only the two values "root", "parterns". 
        The "root" key references a particular amino acid, "partners" references 
        all its interactors.
        When amino acid are part of the same structure the "root" partner resID 
        is always smaller than the ones of its partners.
        When part of different structure, the chain ID of the root amino acid is 
        always the same.
        This guarantees that a particular amino acid pair is mentioned only once.

        eg: The amino acid A41 forms a single contact with amino acid A36
        { "root":{"resID" : "41 ", "chainID":"A"},
        "partners":[{"resID" : "36 ", "chainID":"A"}]
        }
    """
    
    # Thread worker reference
    wThread = None
    # Setting default thread positional arguments
    threadArgs = []
    # Setting default thread named arguments
    threadKwargs = {\
        "encode" : kwargs['encode'] if "encode" in kwargs else False,\
        "atomic" : kwargs['atomic'] if "atomic" in kwargs else False,\
        "deserialize" : kwargs['deserialize']\
                        if "deserialize" in kwargs\
                        else False\
    } 
    
    if isTransformAsFile(kwargs):
    # Create Thread parameters for list of transformation
        wThread                   = lzThread
        threadNum, threadArgs, _threadKwargs = lzGenerateArgsFromFile(\
                                        kwargs['transformation'],\
                                        kwargs['pdbA'],\
                                        kwargs['pdbB'],\
                                        threadNum)
        threadKwargs.update(_threadKwargs)
    # Create Thread parameters for list of PDB structures objects
    elif isExplicitAsFile(kwargs):
        wThread    = lcThread
         threadNum, threadArgs = lcGenerateArgsFromFile(\
                                        kwargs['structList'],\
                                        threadNum)
    elif isTransformAsVector(kwargs):
        wThread                   = lzThread
        threadNum, threadArgs, _threadKwargs = lzGenerateArgsFromVector(\
                                        kwargs['pdbA'],\
                                        kwargs['pdbB'],\
                                        kwargs['euler'],\
                                        kwargs['translation'],\
                                        kwargs['offsetRec'],\
                                        kwargs['offsetLig'],\
                                        threadNum)
    elif isExplicitAsVector(kwargs):
        wThread    = lcThread
        args =       [ kwargs['pdbAtomList'] ]\
                if "pdbAtomList" in kwargs\
                else [ kwargs['pdbAtomListREC'], kwargs['pdbAtomListREC'] ]
        threadNum, threadArgs, _threadKwargs = lcGenerateArgsFromVector(\
                                    *args, threadNum=threadNum)
    else:
        raise TypeError("Unspecified multithreading type")
    
    mStart     = time.time()
    output     = [ None for i in range(threadNum) ]
    threadPool = [  threading.Thread(\
                              args = tuple( [ threadArgs[i], i, output ] )\
                          , kwargs = threadKwargs \
                          , target = wThread) \
                    for i in range(threadNum)\
                 ]

    for th in threadPool:
        th.start()

    for th in threadPool:
        th.join()

    print(f"{threadNum} threads finished in { time.time() - mStart }")        
    
    return output

def lcThread(*args, deserialize=False, **kwargs): 
    tStart = time.time()
    tArgs, tNum, results = args
    assert len(tArgs) == 2
  

    #results[tNum] = []
    print(f"Starting lcThread {tNum}")
    print( f"Shapes of ligand/receptor arrays {len(tArgs[0])}/{len(tArgs[1])}" )
    
    #deserialize = kwargs.pop("deserialize")
    _ = ccmap.lcmap(*tArgs, **kwargs)
    if deserialize:
        results[tNum] = json.loads(_)
    else:
        results[tNum] = _

    print(f"End of lcThread {tNum} in { time.time() - tStart }")          
    return

def lzThread(*args, deserialize=False, **kwargs):
    tStart = time.time()
    tArgs, tNum, results = args
    print(f"Starting lzThread {tNum}")
    #deserialize = kwargs.pop("deserialize")
    _ = ccmap.lzmap(*tArgs, **kwargs)
    if deserialize:
        results[tNum] = json.loads(_)
    else:
        results[tNum] = _
        
    print(f"End of lzThread {tNum} in { time.time() - tStart }")          
    return