# Functions for getting results from output files.

import numpy as np
import json
import dendropy
from copy import deepcopy
import re
from itertools import combinations

from dendropy.calculate import treecompare


def get_value_after_print_line(oneFile, printStr, preStr=None, emptyLineBound=5):
    emptyLine = 0
    f = open(oneFile)
    oneLine = f.readline().rstrip('\n')
    if preStr is not None:
        while (oneLine != preStr) and oneLine:
            oneLine = f.readline().rstrip('\n')
    while (oneLine != printStr):
        oneLine = f.readline().rstrip('\n')
        if oneLine == printStr:
            res = float(f.readline().rstrip('\n'))
            return res
            break
        elif not oneLine:
            emptyLine += 1
            if emptyLine >= emptyLineBound:
                print 'Error String not found!'
                return None
                break
        else:
            emptyLine = 0


def get_value_after_print_line_files(files, printStr, preStr=None):
    res = []
    for oneFile in files:
        res.append(get_value_after_print_line(oneFile, printStr, preStr))
    return res


def switch_out_key_with_in_key_for_dict_in_dict(oneDict):
    res = {}
    for kOut, vOut in oneDict.iteritems():
        for kIn, vIn in vOut.iteritems():
            if kIn in res:
                res[kIn] = {}
            res[kIn][kOut] = vIn
    return res


def len_segment_cut_into_n_subset(lenSegs, lenSubs):
    res = []
    while len(lenSubs) > 0:
        lenSegsCumSum = np.array(lenSegs).cumsum()
        lenSubsCumSum = np.array(lenSubs).cumsum()
        cutPointIndex = [(lenSegsCumSum > cutValue).argmax() for cutValue in lenSubsCumSum]
        lenSegsStartIndex = 0
        cutPoint = cutPointIndex[0]
        lenSegsEndIndex = cutPoint
        tem = lenSegs[lenSegsStartIndex:lenSegsEndIndex]
        if cutPoint == 0:
            reminder = lenSubs[0]
        else:
            reminder = lenSubsCumSum[0] - lenSegsCumSum[cutPoint-1]
        if reminder != 0:
            tem.append(reminder)
        lenSegs = lenSegs[lenSegsEndIndex:]
        lenSegs[0] = lenSegs[0] - reminder
        lenSubs = lenSubs[1:]
        res.append(tem)
    return res


def seg_rate_dict_cut_into_n_subset(segRateDict, lenSegs, lenSubs):
    res = []
    segRateDictValues = segRateDict.values()
    while len(lenSubs) > 0:
        segRateDictNew = {}
        lenSegsCumSum = np.array(lenSegs).cumsum()
        lenSubsCumSum = np.array(lenSubs).cumsum()
        cutPointIndex = [(lenSegsCumSum > cutValue).argmax() for cutValue in lenSubsCumSum]
        lenSegsStartIndex = 0
        cutPoint = cutPointIndex[0]
        lenSegsEndIndex = cutPoint
        tem = lenSegs[lenSegsStartIndex:lenSegsEndIndex]
        if cutPoint == 0:
            reminder = lenSubs[0]
            segRateDictNew[0] = segRateDictValues[0]
        else:
            reminder = lenSubsCumSum[0] - lenSegsCumSum[cutPoint-1]
            for segId in xrange(cutPoint):
                segRateDictNew[segId] = segRateDictValues[segId]
        if reminder != 0:
            tem.append(reminder)
            segRateDictNew[cutPoint] = segRateDictValues[cutPoint]
        segRateDictValues = segRateDictValues[cutPoint:]
        lenSegs = lenSegs[lenSegsEndIndex:]
        lenSegs[0] = lenSegs[0] - reminder
        lenSubs = lenSubs[1:]
        res.append(segRateDictNew)
    return res


def read_all_dicts(dictFiles):
    res = []
    for dictFile in dictFiles:
        with open(dictFile) as json_data:
            res.append(json.load(json_data))
    return res


def dist_tree_all(treeFiles, treeTrueFile):
    treeTrue = dendropy.Tree.get_from_path(treeTrueFile, schema='newick')
    treeTreeTotalLength = treeTrue.length()
    treeTrueScaled = deepcopy(treeTrue)
    treeTrueScaled.scale_edges(1./treeTreeTotalLength)
    distRf = []
    distRfScaled = []
    distSym = []
    for treeFile in treeFiles:
        tree = dendropy.Tree.get_from_path(treeFile, schema='newick')
        distSym.append(treecompare.symmetric_difference(treeTrue, tree))
        distRf.append(treecompare.weighted_robinson_foulds_distance(treeTrue, tree))
        tree.scale_edges(1. / tree.length())
        distRfScaled.append(treecompare.weighted_robinson_foulds_distance(treeTrueScaled, tree))
    return distRf, distRfScaled, distSym


def raw_dicts(dicts):
    oneDict = dicts[0]
    res = {}
    for ka in oneDict.keys():
        res[ka] = {}
        oneDictOneLine = oneDict[ka]
        for kb in oneDictOneLine.keys():
            vec = np.array([dicti[ka][kb] for dicti in dicts])
            res[ka][kb] = vec
    return res


def mean_dict(rawDict):
    res = {}
    for ka, va in rawDict.iteritems():
        res[ka] = {}
        for kb, vb in va.iteritems():
            res[ka][kb] = np.mean(rawDict[ka][kb])
    return res


def median_dict(rawDict):
    res = {}
    for ka, va in rawDict.iteritems():
        res[ka] = {}
        for kb, vb in va.iteritems():
            res[ka][kb] = np.median(rawDict[ka][kb])
    return res


def ctmc0_multi_align(multiAlign, cList, missingSym='?'):
    """
    prepare multiple alignment for CTMC0 model, i.e, no substitution considered
    replace all non-empyt characters by a missing symbol
    input:
        multiAlign: dict, multiple alignment
        missingSym: a symbol to replace all non-empty characters, by default ?
    output:
        multiAlignCTMC0: dict
    """
    multiAlignCTMC0 = {}
    cListJoined = ''.join(cList)
    cListJoinedRe = '[' + cListJoined + ']'
    for k, v in multiAlign.iteritems():
        vNew = re.sub(cListJoinedRe, missingSym, v)
        multiAlignCTMC0[k] = vNew
    return multiAlignCTMC0


def dist_among_trees(treeDict):
    """
    distance matrix of Robinson Foulds difference between every pair of trees
    """
    res = {}
    for treeName1 in treeDict.keys():
        tree1 = treeDict[treeName1]
        res[treeName1] = {}
        for treeName2 in treeDict.keys():
            tree2 = treeDict[treeName2]
            res[treeName1][treeName2] = treecompare.weighted_robinson_foulds_distance(tree1, tree2)
    return res


def dist_among_trees_sym(treeDict):
    """
    distance matrix of symmetric difference between every pair of trees
    """
    res = {}
    for treeName1 in treeDict.keys():
        tree1 = treeDict[treeName1]
        res[treeName1] = {}
        for treeName2 in treeDict.keys():
            tree2 = treeDict[treeName2]
            res[treeName1][treeName2] = treecompare.symmetric_difference(tree1, tree2)
    return res


def all_dist_among_trees(treeDict):
    """
    distance matrix of Robinson Foulds difference between every pair of trees
    """
    res = []
    keys = treeDict.keys()
    comb = combinations(keys, 2)
    for treeName1, treeName2 in comb:
        tree1 = treeDict[treeName1]
        tree2 = treeDict[treeName2]
        res.append(treecompare.weighted_robinson_foulds_distance(tree1, tree2))
    return res


def all_dist_among_trees_sym(treeDict):
    """
    distance matrix of Robinson Foulds difference between every pair of trees
    """
    res = []
    keys = treeDict.keys()
    comb = combinations(keys, 2)
    for treeName1, treeName2 in comb:
        tree1 = deepcopy(treeDict[treeName1])
        tree2 = deepcopy(treeDict[treeName2])
        res.append(treecompare.symmetric_difference(tree1, tree2))
    return res


def output_dist_dict(distDict, outFile):
    """
    distance matrix of Robinson Foulds difference between every pair of trees
    """
    f = open(outFile, 'w')
    keys = distDict.keys()
    head = '\t' + '\t'.join(keys) + '\n'
    f.write(head)
    for oneKey in keys:
        distDictOneLine = distDict[oneKey]
        distVec = distDictOneLine.values()
        oneLine = oneKey + '\t' + '\t'.join([value.__str__() for value in distVec]) + '\n'
        f.write(oneLine)
    f.close()
