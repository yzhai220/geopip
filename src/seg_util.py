# Utility functions for segmentation.

from copy import deepcopy
import dendropy
import numpy as np


def subtree_retain_taxa_with_labels(treeFull, labels):
    """
    Get a subtree with taxa lables.
    """
    treeFullCopy = deepcopy(treeFull)
    treeFullCopy.retain_taxa_with_labels(labels)
    treeFullCopy.suppress_unifurcations()
    treeStr = treeFullCopy.as_string(schema='newick')
    index1 = treeStr.rfind(':')
    index2 = treeStr.rfind(')')
    if index1 > index2:
        treeStrSub = treeStr[(treeStr.find(' ')+1):index1]
    else:
        treeStrSub = treeStr[(treeStr.find(' ')+1):treeStr.rfind(';')]
    treeStrSub = treeStrSub + ';'
    treeSub = dendropy.Tree.get_from_string(treeStrSub, schema='newick')
    return treeSub


def msa_retain_taxa_with_labels(multiAlign, labels, lenSegs):
    """
    Get a MSA from a subtree with taxa labels.
    """
    nLabels = len(labels)
    multiAlignSub = {}
    for label in labels:
        multiAlignSub[label] = multiAlign[label]
    msaListSub = zip(*multiAlignSub.values())
    keysAll = multiAlignSub.keys()
    start = 0
    lenSegsSub = []
    msaListSubEmptyRemoved = []
    for lenSeg in lenSegs:
        end = start + lenSeg
        msaListSubSeg = msaListSub[start:end]
        emptyElement = '- ' * nLabels
        emptyElement = tuple(emptyElement.split())
        nEmpty = msaListSubSeg.count(emptyElement)
        lenSegsSub.append(lenSeg-nEmpty)
        msaListSubSegEmptyRemoved = [msa for msa in msaListSubSeg if msa != emptyElement]
        msaListSubEmptyRemoved += msaListSubSegEmptyRemoved
        start = end
    msaListSub = zip(*msaListSubEmptyRemoved)
    msaListSub = [''.join(value) for value in msaListSub]
    multiAlignSub = dict(zip(keysAll, msaListSub))
    return multiAlignSub, lenSegsSub


def msa_retain_taxa_with_labels_and_track_col_num(multiAlign, labels, lenSegs):
    """
    Track column number.
    """
    nLabels = len(labels)
    multiAlignSub = {}
    # multiAlignAll = {}
    for label in labels:
        multiAlignSub[label] = multiAlign[label]
    msaListSub = zip(*multiAlignSub.values())
    # msaListAll = zip(*multiAlign.values())
    keysAll = multiAlignSub.keys()
    start = 0
    lenSegsSub = []
    msaListSubEmptyRemoved = []
    # msaListAllEmptyRemoved = []
    msaListIndexNonEmtpy = []
    for lenSeg in lenSegs:
        end = start + lenSeg
        msaListSubSeg = msaListSub[start:end]
        # msaListAllSeg = msaListAll[start:end]
        emptyElement = '- ' * nLabels
        emptyElement = tuple(emptyElement.split())
        nEmpty = msaListSubSeg.count(emptyElement)
        lenSegsSub.append(lenSeg-nEmpty)
        msaListSubSegEmptyRemoved = []
        # msaListAllSegEmptyRemoved = []
        # msaListSubSegEmptyRemoved = [msa for msa in msaListSubSeg if msa!=emptyElement]
        for i in xrange(len(msaListSubSeg)):
            msaSub = msaListSubSeg[i]
            # msaAll = msaListAllSeg[i]
            if msaSub != emptyElement:
                msaListSubSegEmptyRemoved.append(msaSub)
                # msaListAllSegEmptyRemoved.append(msaAll)
                msaListIndexNonEmtpy.append(i+start)
        msaListSubEmptyRemoved += msaListSubSegEmptyRemoved
        # msaListAllEmptyRemoved += msaListAllSegEmptyRemoved
        start = end
    msaListSub = zip(*msaListSubEmptyRemoved)
    msaListSub = [''.join(value) for value in msaListSub]
    multiAlignSub = dict(zip(keysAll, msaListSub))
    # msaListAll = zip(*msaListAllEmptyRemoved)
    # msaListAll = [''.join(value) for value in msaListAll]
    # multiAlignAll = dict(zip(multiAlign.keys(), msaListAll))
    return multiAlignSub, lenSegsSub, msaListIndexNonEmtpy


def indel_rate_for_all_loci(segRateDict, ratesList, lenSegs):
    """
    Get indel rate for all loci.
    """
    segIds = segRateDict.keys()
    segIds.sort()
    i = 0
    res = []
    for segId in segIds:
        segRate = segRateDict[segId]
        segRateIndex = ratesList.index(segRate)
        res += [segRateIndex] * lenSegs[i]
        i += 1
    return res


def msa_subset_selection(msaColIndex, msaColIndexSmall):
    """
    return a list of index for loci of elements in msaColIndexSmall with respect to msaColIndex
    """
    res = [msaColIndex.index(msaColIndexOne) for msaColIndexOne in msaColIndexSmall]
    return res


def prop_of_col_with_wrong_rates(indelRatesAll, msaColIndexSmall, msaColIndex, lenSegsSubEst, rateSegsSubEst):
    """
    Calcluate proportion of columns with wrong indel rates.
    """
    indelRatesSubTrue = []
    for index in msaColIndexSmall:
        indelRatesSubTrue += [indelRatesAll[index]]
    indelRatesAllEst = []
    for i in xrange(len(lenSegsSubEst)):
        indelRatesAllEst += [rateSegsSubEst[i]] * lenSegsSubEst[i]
    indelRatesSubEst = []
    for i in xrange(len(indelRatesSubTrue)):
        indelRatesSubEst += [indelRatesAllEst[msaColIndex.index(msaColIndexSmall[i])]]
    ratio = float(np.sum(np.array(indelRatesSubTrue) != np.array(indelRatesSubEst))) / len(indelRatesSubEst)
    return ratio
