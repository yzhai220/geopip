# Geo PIP for segments simulation.

import numpy as np
import random as rd
from bisect import bisect
from copy import deepcopy
import itertools
import dendropy

from pip_util import pi_from_qmat


dnaList = ['A', 'C', 'G', 'T']


def sim_segs_initial(p, ratesList, piProbSeg, piProb, cList, fixSegNumber=False):
    """
    initialize segments
    input:
        p: parameter for geometric distribution
        ratesList: a list of possible rates
        piProbSeg: stationary distribution for all possible rates
        piProb: stationary distributio of characters
        cList: list of characters from a finite alphabet
    output:
        res: a dictionary of results simulated, including iRate, dRate, dnaSeq and locDict
    """
    res = {}
    if fixSegNumber is False:
        nSeg = np.random.geometric(p)
    else:
        nSeg = fixSegNumber
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
        res[locSeg] = {'iRate': iRate, 'dRate': dRate, 'dnaSeq': dnaSeq, 'locDict': locDict}
        i += 1
    return res


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


def segs_change(segs, t, qMat, piProb, cList=dnaList):
    """
    generate one DNA sequence from segs, after time t
    input:
        segs: from sim_segs_initial
        t: evolutionary distance
        qMat, piProb: rate matrix and stationary distribution
    output:
        segsNew: new segment, a dictionary
    """
    cumProb = np.array(piProb).cumsum()
    segsNew = deepcopy(segs)
    for segLoc, segValue in segsNew.iteritems():
        segValue['dnaSeq'], segValue['locDict'] = in_seg_change_pip(segValue['dnaSeq'], segValue['locDict'], piProb, cumProb, segValue['iRate'], segValue['dRate'], qMat, t, cList)
    return segsNew


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
        subRateVec = np.array([-qMat[cList.index(ch), cList.index(ch)] for ch in word])
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
        elif tMinIndex == 1:
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
            # print ch
            chInd = cList.index(ch)    # index of the old character
            chSubRate = qMat[chInd, ].copy()   # transition rate
            chSubRate = chSubRate.clip(min=0)
            chSubRateScaled = chSubRate / chSubRate.sum()   # scaled to 1
            chNewVec = np.random.multinomial(1, chSubRateScaled)
            chNewInd = chNewVec.argmax()    # find new character index
            chNew = cList[chNewInd]    # find new character
            # print chnew
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


def sim_tree(tree, p, ratesList, piProbSeg, piProb, qMat, cList, fixSegNumber=False):
    """
    Generate sequences based on a tree using the GeoPIP model.
    input:
        tree: a tree, of dendropy.Tree class
        p: parameter of geometric distribution
        ratesList: a list of all rates
        piProb, qMat, cList: ...
    output:
        tree: updated tree, with value at each node
    """
    piProb = pi_from_qmat(qMat)
    for node in tree.preorder_node_iter():
        if node.parent_node is None:
            value = sim_segs_initial(p, ratesList, piProbSeg, piProb, cList, fixSegNumber)
            node.value = value
        else:
            value = segs_change(node.parent_node.value, node.edge.length, qMat, piProb, cList)
            node.value = value
    return tree


def leaf_data_from_sim_tree(tree):
    """
    summarizing results from simulated tree
    output:
    """
    seqs = {}
    alignInSeg = {}
    alignSeg = {}
    segRateDict = {}
    # notAlignSeqs = {}
    ts = {}
    treeTem = deepcopy(tree)
    pdm = dendropy.treecalc.PatristicDistanceMatrix(treeTem)
    leafNodes = tree.leaf_nodes()
    lanNames = [leafNode.taxon.label for leafNode in leafNodes]
    lanNames.sort()
    lanPairsIter = itertools.combinations(lanNames, 2)
    segRateDict = get_seg_rate_dict(lanNames, tree)
    segIds = segRateDict.keys()
    segIds.sort()
    nSegs = len(segIds)
    segIdsNew = range(nSegs)
    segIdsTrans = dict(zip(segIds, segIdsNew))
    for segId in segIds:
        segIdNew = segIdsTrans[segId]
        segRateDict[segIdNew] = segRateDict.pop(segId)
    # get simple sequences without alignment
    for lanName in lanNames:
        node = tree.find_node_with_taxon_label(lanName)
        seqs[lanName] = node.value
        # notAlignSeqs[lanName] = get_dna_seq(node.value)
    # get multiple alignment
    segIdLocIdDict = get_all_locid_within_all_segs(seqs, segIds)
    multiAlignSeqs = get_dna_multi_align(segIds, segIdLocIdDict, seqs, lanNames)
    for lanName, values in seqs.iteritems():
        segIdsInLan = values.keys()
        for segId in segIdsInLan:
            segIdNew = segIdsTrans[segId]
            seqs[lanName][segIdNew] = seqs[lanName][segId]
            del seqs[lanName][segId]
    # get pairwise alignemtn
    for lanPair in lanPairsIter:
        lanName1 = lanPair[0]
        lanName2 = lanPair[1]
        node1 = tree.find_node_with_taxon_label(lanName1)
        node2 = tree.find_node_with_taxon_label(lanName2)
        segs1 = node1.value
        segs2 = node2.value
        alignInSeg[lanPair] = get_dna_align(segs1, segs2)
        alignSeg[lanPair] = get_seg_align(segs1, segs2)
        ts[lanPair] = pdm(node1.taxon, node2.taxon)
        # very small compuatational error in pdm method
    # get segment lengths
    seq0 = deepcopy(multiAlignSeqs[multiAlignSeqs.keys()[0]])
    lenSegs = []
    while seq0 != '':
        next = seq0.index('*')
        if next != 0:
            lenSegs.append(next)
            seq0 = seq0[(next+1):]
        else:
            seq0 = seq0[1:]
    return seqs, alignInSeg, alignSeg, segRateDict, ts, multiAlignSeqs, lenSegs


def get_seg_rate_dict(leafNames, tree):
    """
    get all IDs of segments from leafs
    """
    segRateDict = {}
    for leafName in leafNames:
        node = tree.find_node_with_taxon_label(leafName)
        nodeValue = node.value
        for k, v in nodeValue.iteritems():
            if k not in segRateDict:
                segRateDict[k] = (v['iRate'], v['dRate'])
    return segRateDict


def get_dna_seq(segs):
    """
    summarize DNA sequence data of one species
    """
    dnaSeq = '*'   # symble * separate DNA sequences from different segment
    segLocList = segs.keys()
    segLocList.sort()
    for segLoc in segLocList:
        dnaSeq = dnaSeq + segs[segLoc]['dnaSeq'] + '*'
    return dnaSeq


def get_seg_align(segs1, segs2):
    """
    get true alignment of segments of two taxa
    """
    alignSeg = []
    segs1.keys()
    seg1LocList = segs1.keys()
    seg1LocList.sort()
    seg2LocList = segs2.keys()
    seg2LocList.sort()
    segLocListUnion = list(set(seg1LocList + seg2LocList))
    segLocListUnion.sort()
    for segLoc in segLocListUnion:
        boolTem1 = (segLoc in segs1)
        boolTem2 = (segLoc in segs2)
        if boolTem1 and boolTem2:
            rate1 = (segs1[segLoc]['iRate'], segs1[segLoc]['dRate'])
            rate2 = (segs2[segLoc]['iRate'], segs2[segLoc]['dRate'])
            alignSeg += [(rate1, rate2)]
        elif boolTem1 and not boolTem2:
            rate1 = (segs1[segLoc]['iRate'], segs1[segLoc]['dRate'])
            rate2 = '-'
            alignSeg += [(rate1, rate2)]
        else:
            rate1 = '-'
            rate2 = (segs2[segLoc]['iRate'], segs2[segLoc]['dRate'])
            alignSeg += [(rate1, rate2)]
    return alignSeg


def get_dna_align(segs1, segs2):
    """
    get true alignment of two DNA sequences from simulation for all segments
    """
    align = {}
    seg1LocList = segs1.keys()
    seg1LocList.sort()
    seg2LocList = segs2.keys()
    seg2LocList.sort()
    segLocListUnion = list(set(seg1LocList + seg2LocList))
    segLocListUnion.sort()
    for segLoc in segLocListUnion:
        boolTem1 = (segLoc in segs1)
        boolTem2 = (segLoc in segs2)
        if boolTem1 and boolTem2:
            locDict1 = segs1[segLoc]['locDict']
            locDict2 = segs2[segLoc]['locDict']
            align[segLoc] = align_from_loc_one_word(locDict1, locDict2)
        elif boolTem1 and not boolTem2:
            locDict1 = segs1[segLoc]['locDict']
            locDict2 = {}
            align[segLoc] = align_from_loc_one_word(locDict1, locDict2)
        else:
            locDict1 = {}
            locDict2 = segs2[segLoc]['locDict']
            align[segLoc] = align_from_loc_one_word(locDict1, locDict2)
    return align


def align_from_loc_one_word(locDict1, locDict2):
    """
    reconstruct true alignment based on 2 locDict by matching keys
    only one word is considered, with only values, no word id
    """
    res = []
    combKeys = list(set(locDict1.keys() + locDict2.keys()))
    combKeysSorted = sorted(combKeys)
    for key in combKeysSorted:
        if key in locDict1:
            va = locDict1[key]
        else:
            va = '-'
        if key in locDict2:
            vb = locDict2[key]
        else:
            vb = '-'
        res.append((va, vb))
    return res


def get_all_locid_within_one_seg(seqs, segId):
    locIdAllInSegId = []
    for seqId, seqValue in seqs.iteritems():
        if segId in seqValue:
            locIdAllInSegId += seqValue[segId]['locDict'].keys()
    locIdAllInSegId = list(set(locIdAllInSegId))
    locIdAllInSegId.sort()
    return locIdAllInSegId


def get_all_locid_within_all_segs(seqs, segIds):
    segIdLocIdDict = {}
    for segId in segIds:
        segIdLocIdDict[segId] = get_all_locid_within_one_seg(seqs, segId)
    return segIdLocIdDict


def get_dna_multi_align(segIds, segIdLocIdDict, seqs, lanNames):
    """
    get multiple alignment from simulated data
    """
    multiAlignSeqs = {}
    segIds.sort()
    for lan in lanNames:
        lanValue = seqs[lan]
        seqAllSegs = ''
        for segId in segIds:
            seqInSeg = '*'
            if segId in lanValue:
                lanInSegIdLocDict = lanValue[segId]['locDict']
                locIdAllSorted = segIdLocIdDict[segId]
                for locId in locIdAllSorted:
                    if locId in lanInSegIdLocDict:
                        seqInSeg += lanInSegIdLocDict[locId]
                    else:
                        seqInSeg += '-'
            else:
                seqInSeg += '-' * len(segIdLocIdDict[segId])
            seqAllSegs += seqInSeg
        multiAlignSeqs[lan] = seqAllSegs + '*'
    return multiAlignSeqs
