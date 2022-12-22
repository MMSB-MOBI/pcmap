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


# Alternative PoolExecutor

def run_over_np_frame(frame_as_np_arr):
    pass