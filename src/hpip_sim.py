# The hPIP for segments simulation.

import numpy as np
import random as rd
from bisect import bisect
from copy import deepcopy
from pip_util import pi_from_qmat

dnaList = ['A', 'C', 'G', 'T']


def sim_segs_initial(iRateSeg, dRateSeg, piProbSeg, ratesList, piProb, cList=dnaList, fixSegNumber=False):
    """
    initialize segments
    """
    if fixSegNumber == False:
        nSeg = np.random.poisson(iRateSeg / dRateSeg)
    else:
        nSeg = fixSegNumber
    res = {}
    locSegVec = np.random.uniform(size=nSeg)
    locSegVec.sort()
    cumProbSeg = np.array(piProbSeg).cumsum()
    segTypeSeed = np.random.uniform(size=nSeg)
    i = 0
    for locSeg in locSegVec:
        segType = segTypeSeed[i]
        rate = ratesList[bisect(cumProbSeg, segType)]
        iRate = rate[0]
        dRate = rate[1]
        dnaSeq, locDict = sim_one_seg_initial(iRate, dRate, piProb, cList)
        res[locSeg] = {'iRate':iRate, 'dRate':dRate, 'dnaSeq':dnaSeq, 'locDict':locDict}
        i += 1
    return res


# same as the one in sim_geopip_2seq.py
def sim_one_seg_initial(iRate, dRate, piProb, cList=dnaList):
    """
    simulate sequence data in one segment based on PIP
    input:
        iRate, dRate: insertion and deletion rate
        piProb: stationary distribution of characters
        cList: list of characters
    output:
        res: a string of characters simulated
        locDict: a dictionary of location->character
    """
    segLen = np.random.poisson(iRate / dRate)
    locDict = {}
    cumProb = np.array(piProb).cumsum()
    locVec = np.random.uniform(size=segLen)
    locVec.sort()
    chSeeds = np.random.uniform(size=segLen)
    chAllList = [cList[bisect(cumProb, chSeed)] for chSeed in chSeeds]
    res = ''.join(chAllList)   # we can remove res later, not used in inference
    locDict = dict(zip(locVec, chAllList))
    return res, locDict


def segs_change(segs, piProbSeg, iRateSeg, dRateSeg, qMatSeg, t, ratesList, qMat, piProb, cList=dnaList):
    """generate one DNA sequence"""
    cumProbSeg = np.array(piProbSeg).cumsum()
    cumProb = np.array(piProb).cumsum()
    segsNew = out_seg_change_pip(segs, piProbSeg, cumProbSeg, iRateSeg, dRateSeg, qMatSeg, t, ratesList, piProb, cList)
    for segLoc, segValue in segsNew.iteritems():
        segValue['dnaSeq'], segValue['locDict'] = in_seg_change_pip(segValue['dnaSeq'], segValue['locDict'], piProb, cumProb, segValue['iRate'], segValue['dRate'], qMat, t, cList)
    return segsNew


def out_seg_change_pip(segs, piProbSeg, cumProbSeg, iRateSeg, dRateSeg, qMatSeg, t, ratesList, piProb, cList=dnaList):
    """generate one DNA sequence in one segment based on PIP"""
    stop = False
    segsNew = deepcopy(segs)
    tNew = deepcopy(t)
    while not stop:
        segsNew, tNew, stop = out_seg_change_pip_iter(segsNew, piProbSeg, cumProbSeg, iRateSeg, dRateSeg, qMatSeg, tNew, ratesList, piProb, cList)
        # print segsNew.keys()
    return segsNew


def out_seg_change_pip_iter(segs, piProbSeg, cumProbSeg, iRateSeg, dRateSeg, qMatSeg, t, ratesList, piProb, cList=dnaList):
    """
    one iteration only
    generate one new DNA sequence based on PIP
    return: wordNew, stop (indicate whether time is bigger than t), t, cType, locDict
    """
    nSeg = len(segs)
    tIns = np.random.exponential(scale=1./iRateSeg)
    # rateAllCurrent = [(segValue['iRate'], segValue['dRate']) for segLoc, segValue in segs.iteritems()]
    if nSeg == 0:    # if word is empty, only insertion can happen
        tMinIndex = 0
        tMin = tIns
    else:   # if word is not empty
        tDel = np.array(np.random.exponential(scale=1./dRateSeg, size=nSeg))
        # only INSERTION and DELETION of segments
        # set substitution rate to be small, 1.e-10
        subRateVec = np.zeros(nSeg) + 1.e-10
        tSub = np.array(np.random.exponential(scale=1./subRateVec))
        # print tSub
        tDelMin = tDel.min()
        tDelMinIndex = tDel.argmin()
        tSubMin = tSub.min()
        # print tSubIndex
        tVec = np.array([tIns, tDelMin, tSubMin])
        tMinIndex = tVec.argmin()
        tMin = tVec.min()
    if tMin < t:
        if tMinIndex == 0:    # insertion
            segNewLoc = rd.random()    # location
            rateSeed = rd.random()
            rate = ratesList[bisect(cumProbSeg, rateSeed)]
            iRate = rate[0]
            dRate = rate[1]
            dnaSeq, locDict = sim_one_seg_initial(iRate, dRate, piProb, cList)
            segNewValue = {'iRate':iRate, 'dRate':dRate, 'dnaSeq':dnaSeq, 'locDict':locDict}
            segs[segNewLoc] = segNewValue
            cType = 'insertion'
        elif tMinIndex ==1:
            # deletion
            segs.pop(segs.keys()[tDelMinIndex])
            cType = 'deleletion'
        else:
            # substitution
            cType = 'deleletion'
        t = t - tMin
        stop = False
    else:
        stop = True
        cType = 'none'
    # print cType
    return segs, t, stop



def in_seg_change_pip(word, locDict, piProb, cumProb, iRate, dRate, qMat, t, cList):
    """generate one DNA sequence in one segment based on PIP"""
    stop = False
#    print word, t
    while not stop:
        word, stop, t, cType, locDict = in_seg_change_pip_iter(word, locDict, piProb, cumProb, iRate, dRate, qMat, t, cList)
    return word, locDict


def in_seg_change_pip_iter(word, locDict, piProb, cumProb, iRate, dRate, qMat, t, cList):
    """
    one iteration only
    generate one new DNA sequence based on PIP
    return: wordNew, stop (indicate whether time is bigger than t), t, cType, locDict
    """
    wordLen = len(word)
    tIns = np.random.exponential(scale=1./iRate)
    if wordLen == 0:    # if word is empty, only insertion can happen
        tMinIndex = 0
        tMin = tIns
    else:   # if word is not empty
        tDel = np.array(np.random.exponential(scale=1./dRate, size=wordLen))
        subRateVec = np.array([-qMat[cList.index(ch),cList.index(ch)] for ch in word])
        tSub = np.array(np.random.exponential(scale=1./subRateVec))
        # print tSub
        tDelMin = tDel.min()
        tSubMin = tSub.min()
        tDelIndex = tDel.argmin()
        tSubIndex = tSub.argmin()
        # print tSubIndex
        tVec = np.array([tIns, tDelMin, tSubMin])
        tMinIndex = tVec.argmin()
        tMin = tVec.min()
    if tMin < t:
        if tMinIndex == 0:    # insertion
            locVal = rd.random()    # location
            loc = np.sum(np.array(locDict.keys()) < locVal)
            chRan = rd.random()   # character
            chInd = bisect(cumProb, chRan)
            chNew = cList[chInd]    # find new character
            wordNew = word[:loc] + chNew + word[loc:]    # gen new word
            locDict[locVal] = chNew
            cType = 'ins'
        elif tMinIndex ==1:
            # deletion
            loc = tDelIndex    # location
            wordNew = word[:loc] + word[(loc+1):]    # gen new word
            locVal = sorted(locDict.keys())[loc]
            locDict.pop(locVal)
            cType = 'del'
        else:
            # substitution
            loc = tSubIndex    # location
            ch = word[loc]    # character to be replaced
#            print ch
            chInd = cList.index(ch)    # index of the old character
            chSubRate = qMat[chInd, ].copy()   # transition rate
            chSubRate = chSubRate.clip(min=0)
            chSubRateScaled = chSubRate / chSubRate.sum()   # scaled to 1
            chNewVec = np.random.multinomial(1, chSubRateScaled)
            chNewInd = chNewVec.argmax()    # find new character index
            chNew = cList[chNewInd]    # find new character
#            print chnew
            wordNew = word[:loc] + chNew + word[(loc+1):]    # gen new word
            locVal = sorted(locDict.keys())[loc]
            locDict[locVal] = chNew
            cType = 'sub'
        t = t - tMin
        stop = False
    else:
        wordNew = word
        stop = True
        cType = 'none'
        loc = None
    # print cType
    return wordNew, stop, t, cType, locDict


def sim_tree_hpip(tree, iRateSeg, dRateSeg, piProbSeg, qMatSeg, ratesList, piProb, qMat, cList, fixSegNumber=False):
    """
    generate sequences based on a tree using hPIP model.
    """
    piProb = pi_from_qmat(qMat)
    for node in tree.preorder_node_iter():
        if node.parent_node is None:
            value = sim_segs_initial(iRateSeg, dRateSeg, piProbSeg, ratesList, piProb, cList, fixSegNumber)
            node.value = value
        else:
            value = segs_change(node.parent_node.value, piProbSeg, iRateSeg, dRateSeg, qMatSeg, node.edge.length, ratesList, qMat, piProb, cList)
            node.value = value
    return tree
