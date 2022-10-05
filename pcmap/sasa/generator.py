from ..generators import splitInterval


#threadNum, threadArgs, _threadKwargs = sasaManyGenerateArgsFromFileList(
# input iterator must have a length ...
def sasaManyGenerateArgsFromVectors( pdbDictorizedList, baseArgs, baseKwargs, threadNum):
    n_struct = len(pdbDictorizedList)
    threadNum = min(threadNum, n_struct)

    threadArgs = []
    threadKwargs =  []
    for x,y in splitInterval(n_struct, threadNum):
        threadArgs.append( ( pdbDictorizedList[x:y], *baseArgs ) )
        threadKwargs.append(baseKwargs)
    return threadArgs, threadKwargs