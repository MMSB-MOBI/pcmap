def splitInterval(iLen, nElem):
    """Yield successive interval boundaries spanning [0, iLen]"""
    assert(nElem <= iLen)
    nWidth = int(iLen/nElem)
    
    for i in range(nElem):
        top = (i+1) * nWidth
        if i == (nElem - 1): #and iLen%nElem != 0:
            top +=  iLen%nElem
        yield (i * nWidth, top)

def lcGenerateArgs(data, oligomerState, threadNum):
    
    threadNum = min(threadNum, len(data))
    threadArgs = []
    for x,y in splitInterval(len(data), threadNum):
        threadArgs.append(( [], [] )) \
            if oligomerState == "dimer" \
            else threadArgs.append(( [] ))
        for i in range(y - x):
            threadArgs[-1][0].append(data[i][0])
            if strType == "dimer":
                threadArgs[-1][1].append(data[i][1])
    return threadNum, threadArgs

def lcGenerateArgsFromFile(fStructList, threadNum):
    """Generate thread constructor parameters for the lc algorithm
    
    :param fStructList: filepath to the list of structures
    :type file: string
    :param threadNum: the number of threads
    :type threadNum: int
    :return: threadArgs
    :rtype: list
    """
     
    strCnt, strList = parseStructFileTabList(fStructList)
    strType = "dimer" if strCnt == 2 else "monomer"
    threadNum, threadArgs = lcGenerateArgs(strList, strType, threadNum)
    return threadNum, threadArgs

def lcGenerateArgsFromVector(*args, threadNum):
    """Generate thread constructor parameters for the lc algorithm
    
    :param fStructList: filepath to the list of structures
    :type file: string
    :param threadNum: the number of threads
    :type threadNum: int
    :return: threadArgs
    :rtype: list
    """

    strType = "dimer" if len(args) == 2 else "monomer"
    strList = args[0] if strType == "monomer" else [ _ for _ in zip(args[:2]) ]
    
    threadNum, threadArgs = lcGenerateArgs(strList, strType, threadNum)
    return threadNum, threadArgs

def lzGenerateArgsFromfile(fTransformation, pdbRec, pdbLig, threadNum):
    """Generate thread constructor parameters for the lz algorithm
    
    :param fTransformation: filepath to the JSON transformation file
    :type file: string
    :param pdbRec: the proteinA aka receptor
    :type pdbRec: pdb object
    :param pdbLig: the proteinB aka receptor
    :type pdbLig: pdb object
    :return: tuple(threadArgs, threadKwargs)
    :rtype: tuple(list, dict)
    """
    eulers, translations,\
    offsetRec, offsetLig = parseTransformationFile(fTransformation)
    
    threadNum, threadArgs, threadKwargs = lzGenerateArgsFromVector(\
                                            pdbRec, pdbLig,\
                                            eulers, translations,\
                                            offsetRec, offsetLig,\
                                            threadNum)

    return threadArgs, threadKwargs

def lzGenerateArgsFromVector(pdbRec, pdbLig,\
                             euler, translation,\
                             offsetRec, offsetLig,\
                             threadNum):
    """Generate thread constructor parameters for the lz algorithm
    
    :param pdbRec: the proteinA aka receptor
    :type pdbRec: pdb object
    :param pdbLig: the proteinB aka receptor
    :type pdbLig: pdb object
    :param euler: euler triplets
    :type euler: vector3[]
    :param translation: translation triplets
    :type translation: vector3[]
    :param offsetRec:  receptor centering vector
    :type offsetRec: vector3
    :param offsetLig:  ligand centering vector
    :type offsetLig: vector3
    :return: tuple(threadArgs, threadKwargs)
    :rtype: tuple(list, dict)
    """
    nTransform = len(eulers)
    threadNum = min(threadNum, nTransform)
    threadKwargs= { "offsetRec" : offsetRec,
                    "offsetLig" : offsetLig
                  }

    threadArgs = []
    pdbDictREC = pdbRec.atomDictorize
    pdbDictLIG = pdbLig.atomDictorize
    for x,y in splitInterval(nTransform, threadNum):
        threadArgs.append( ( pdbDictREC, pdbDictLIG, \
                            eulers[x:y],\
                            translations[x:y]\
                            ) )
    return threadNum, threadArgs, threadKwargs