from .generator import *
import time
import ccmap
import threading

def _setSasaThreadCommonArgs(**_kwargs):
    args = []
    kwargs = {}
    assert 'r' in _kwargs
    kwargs['probe'] = _kwargs['r']
    assert 'atom_radii' in _kwargs 
    args.append(_kwargs['atom_radii'])
    return args, kwargs

def sasaManyThread(*args, **kwargs):
    tStart = time.time()
    tArgs, tNum, results = args
    print(f"Starting sasaManyThread {tNum}")
    results[tNum] = ccmap.sasa_list(*tArgs, **kwargs)
    print(f"End of sasaManyThread {tNum} in { time.time() - tStart }")          
    return
def run( pdbDictorizedAsList, threadNum=8, **kwargs ):
    """ Run threads to process supplied structures 
    
       
    """
  
    # Setting default thread positional arguments
    threadArgs_common, threadKwargs_common = _setSasaThreadCommonArgs(**kwargs)
    threadArgs, threadKwargs = sasaManyGenerateArgsFromVectors(pdbDictorizedAsList, threadArgs_common, threadKwargs_common, threadNum)
    
    mStart     = time.time()
    output     = [ None for i in range(threadNum) ]
    threadPool = [  threading.Thread(\
                            args = tuple( [ threadArgs[i], i, output ] )\
                        , kwargs = threadKwargs[i] \
                        , target = sasaManyThread) \
                    for i in range(threadNum)\
                 ]

    for th in threadPool:
        th.start()

    for th in threadPool:
        th.join()
    
    return output
        ## END OF ..
'''
    outputFmtPipe = lzDeserialiser
    if isTransformAsFile(kwargs):
    # Create Thread parameters for list of transformation extracted from file
        wThread                   = lzThread
        threadNum, threadArgs, _threadKwargs = lzGenerateArgsFromfile(\
                                        kwargs['transformation'],\
                                        kwargs['pdbA'],\
                                        kwargs['pdbB'],\
                                        threadNum)
        threadKwargs.update(_threadKwargs)
    # Create Thread parameters for list of PDB files path extracted from file
    elif isExplicitAsFile(kwargs):
        outputFmtPipe = lcDeserialiser
        wThread               = lcThread
        threadNum, threadArgs = lcGenerateArgsFromFile(\
                                        kwargs['structList'],\
                                        threadNum)
    # Create Thread parameters for list of transformation passed as parameters
    elif isTransformAsVector(kwargs):
        wThread                              = lzThread
        threadNum, threadArgs, _threadKwargs = lzGenerateArgsFromVector(\
                                                kwargs['pdbA'],\
                                                kwargs['pdbB'],\
                                                kwargs['eulers'],\
                                                kwargs['translations'],\
                                                kwargs['offsetRec'],\
                                                kwargs['offsetLig'],\
                                                threadNum)
        threadKwargs.update(_threadKwargs)
    # Create Thread parameters for a single(s) or structure pair(s) passed 
    # as atom vectors
    elif isExplicitAsVector(kwargs):
        outputFmtPipe = lcDeserialiser
        wThread = lcThread
        args    =       [ kwargs['pdbAtomList'] ]\
                if "pdbAtomList" in kwargs\
                else [ kwargs['pdbAtomListREC'], kwargs['pdbAtomListREC'] ]
        threadNum, threadArgs = lcGenerateArgsFromVector(\
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

    print("Output type:", "deserialize" if\
                           threadKwargs['deserialize'] else\
                           "serialized" )
    
    return outputFmtPipe(output)
'''