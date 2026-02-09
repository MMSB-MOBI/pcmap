import json

def lcDeserialiser(threadsOutput):
    outputList = []
    for _ in threadsOutput:
        threadOutput = [ json.loads(x) for x in _ ]
        for ccmapRef in threadOutput:
            outputList.append(ccmapRef['data'])
    return outputList

def lzDeserialiser(threadsOutput):
    outputList = []
    for _ in threadsOutput:
        threadOutput = json.loads(_)
        for ccmap in threadOutput['data']:
            outputList.append(ccmap)
    return outputList