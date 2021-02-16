import sys, threading, json, time, copy, os, ccmap





'''
Perform mutlithreading using supplied transformation
'''
def zThreading(proteinA, proteinB, **kwargs):
    
'''
Perform multithreading using explicit list of coordinates
'''
def lThreading(receptorlist, ligandList):
    pass


def run(threadNum=8):
    mStart  = time.time()
    output  = [ None for i in range(threadNum) ]
    threadPool = [  threading.Thread(args = tuple( [ threadArgs[i], i, output ] )\
                , kwargs = threadKwargs \
                , target=wThread) \
                for i in range(threadNum)]

    for th in threadPool:
        th.start()

    for th in threadPool:
        th.join()

    print(f"{threadNum} threads finished in { time.time() - mStart }")        


def cThread(*args, **kwargs):
    tArgs, tNum, results = args
    dual = len(tArgs) == 2
    results[tNum] = []
    print(f"cThread {tNum}: Starting as dual:{dual}")
    print( f"cThread {tNum}: Shapes of ligand/receptor arrays {len(tArgs[0])}/{len(tArgs[1])}" )
    
    tStart = time.time()
    if dual:
        for rec,lig in zip(tArgs[0], tArgs[1]):
            results[tNum].append( ccmap.cmap(rec, lig, **kwargs) )
    else :
        for mol in tArgs[0]:
            results[tNum].append( ccmap.cmap(mol, **kwargs) )
    
    print(f"cThread {tNum}: Finished in { time.time() - tStart }")          
    return

def lcThread(*args, **kwargs): 
    tStart = time.time()
    tArgs, tNum, results = args
    assert len(tArgs) == 2
    

    results[tNum] = []
    print(f"Starting lcThread {tNum}")
    print( f"Shapes of ligand/receptor arrays {len(tArgs[0])}/{len(tArgs[1])}" )
    
    tStart = time.time()
    results[tNum] = ccmap.lcmap(*tArgs, **kwargs)
    
    print(f"End of lcThread {tNum} in { time.time() - tStart }")          
    return

def zThread(*args, **kwargs):
    tStart = time.time()
    tArgs, tNum, results = args
    pdbRec, pdbLib, eulerList, translationList = tArgs

    print(f"Starting zThread {tNum}, for {len(eulerList)} calls")
    
    results[tNum] = []
    for e,t in zip(eulerList, translationList):
        results[tNum].append( ccmap.zmap(pdbRec, pdbLib, e, t, **kwargs) )
        
    print(f"End of zThread {tNum} in { time.time() - tStart }")          
    return

def lzThread(*args, **kwargs):
    tStart = time.time()
    tArgs, tNum, results = args
    print(f"Starting lzThread {tNum}")
    
    results[tNum] = ccmap.lzmap(*tArgs, **kwargs)
        
    print(f"End of lzThread {tNum} in { time.time() - tStart }")          
    return