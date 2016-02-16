# Functions for alignment operations.

import numpy as np
from itertools import combinations


def pair_align_from_multi_align(multiAlign, lenSegs=None):
    """
    get pairwise alignment dictionary from multiple alignment with seg lengths
    input:
        multiAlign: dict, seqName->str
        lenSegs: list, length of segments, the sum of which is the length of str in multiAlign
    output:
        pairAlignDict: dict, pair->segId->pairwise alignment list
    """
    pairAlignDict = {}
    multiAlignInSeg, segIds = multi_align_sep_by_seg(multiAlign, lenSegs)
    seqNames = multiAlign.keys()
    seqNames.sort()
    pairList = combinations(seqNames, 2)
    for pair in pairList:
        seqName1, seqName2 = pair
        pairAlignDict[pair] = {}
        seq1SegIdStr = multiAlignInSeg[seqName1]
        seq2SegIdStr = multiAlignInSeg[seqName2]
        for segId in segIds:
            str1 = seq1SegIdStr[segId]
            str2 = seq2SegIdStr[segId]
            pairAlign = zip(str1, str2)
            # remove all ('-', '-')
            pairAlignDict[pair][segId] = [align for align in pairAlign if align != ('-', '-')]
    return pairAlignDict


# This function is used for CTMC estimate only.
def pair_align_from_multi_align_subs_only(multiAlign, lenSegs=None):
    """
    get pairwise alignment dictionary from multiple alignment with seg lengths
    input:
        multiAlign: dict, seqName->str
        lenSegs: list, length of segments, the sum of which is the length of str in multiAlign
    output:
        pairAlignDict: dict, pair->segId->pairwise alignment list
    """
    pairAlignDict = {}
    multiAlignInSeg, segIds = multi_align_sep_by_seg(multiAlign, lenSegs)
    seqNames = multiAlign.keys()
    seqNames.sort()
    pairList = combinations(seqNames, 2)
    for pair in pairList:
        seqName1, seqName2 = pair
        pairAlignDict[pair] = {}
        seq1SegIdStr = multiAlignInSeg[seqName1]
        seq2SegIdStr = multiAlignInSeg[seqName2]
        for segId in segIds:
            str1 = seq1SegIdStr[segId]
            str2 = seq2SegIdStr[segId]
            pairAlign = zip(str1, str2)
            # remove all ('-', '-')
            pairAlignDict[pair][segId] = [align for align in pairAlign if align.count('-') == 0]
    return pairAlignDict


def multi_align_sep_by_seg(multiAlign, lenSegs=None):
    """
    multi alignment separated by segs
    input:
        multiAlign: dict, seqName -> string, for all loci
        lenSegs: list, length of segments
    output:
        multiAlignInSeg: dict, seqName -> segId -> string in seg
        segIds: list, Ids of all segments
    """
    if lenSegs is None:
        lenSegs = [len(multiAlign.values()[0])]
    multiAlignInSeg = {}
    nSeg = len(lenSegs)
    segIds = range(nSeg)
    segBounds = np.array(lenSegs).cumsum()
    for seq, multiAlign in multiAlign.iteritems():
        multiAlignInSeg[seq] = {}
        segStart = 0
        for index in xrange(nSeg):
            segEnd = segBounds[index]
            multiAlignInSeg[seq][segIds[index]] = multiAlign[segStart:segEnd]
            segStart = segEnd
    return multiAlignInSeg, segIds


def get_multi_align_in_one_seg(multiAlignInSeg, segId):
    """
    get multiple alignments in one segment with segId
    input:
        multiAlignInSeg: result from multi_align_sep_by_seg, dict, seqName -> segId -> string in seg
        segId
    output:
        multiAlignOneSeg: dict, seqName -> string in seg
    """
    multiAlignOneSeg = {}
    for seqName, seqValue in multiAlignInSeg.iteritems():
        multiAlignOneSeg[seqName] = seqValue[segId]
    return multiAlignOneSeg


def get_multi_align_in_all_seg(multiAlignInSeg, segIds):
    """
    get multiple alignments in all segments with segId
    input:
        multiAlignInSeg: result from multi_align_sep_by_seg, dict, seqName -> segId -> string in seg
        segIds: result from multi_align_sep_by_seg,
    output:
        multiAlignAllSeg: dict, segId -> seqName -> string in seg
    """
    multiAlignAllSeg = {}
    for segId in segIds:
        multiAlignAllSeg[segId] = get_multi_align_in_one_seg(multiAlignInSeg, segId)
    return multiAlignAllSeg


def get_true_seg_rate(segRateDict):
    """
    organize segRateDict so that seg ID is removed, and segs are ordered
    input:
        segRateDict: from simulated data
    output:
        rateSegsList: a list of segment rates, ordered by segments
    """
    keys = segRateDict.keys()
    keys.sort()
    rateSegsList = []
    for key in keys:
        rateSegsList.append(segRateDict[key])
    return rateSegsList


def get_true_seg_rate_index(rateSegsList, ratesList):
    """
    get true segment rate index according to ratesList order
    input:
        rageSegsList: from get_true_seg_rate()
    output:
        rateSegsIndex: a list of segment rates index
    """
    rateSegsIndex = [ratesList.index(rateSeg) for rateSeg in rateSegsList]
    return rateSegsIndex


def multi_align_no_sep(multiAlignSeqs):
    """
    re-organize multiple alignment by locations
    input:
        multiAlignSeqs: from simulated data
    output:
        multiAlignNoSep: separation symbol '*' is removed
    """
    multiAlignNoSep = {}
    for key, value in multiAlignSeqs.iteritems():
        value = value.replace('*', '')
        multiAlignNoSep[key] = value
    return multiAlignNoSep


def multi_align_by_loc(multiAlignNoSep):
    """
    re-organize multiple alignment by locations
    input:
        multiAlignSeqs: from simulated data
    output:
        seqNameAll: name of all sequences, ordered
        multiAlignByLoc: aligned data by location
    """
    seqNameAll = multiAlignNoSep.keys()
    seqNameAll.sort()
    multiAlignByLoc = []
    for seqName in seqNameAll:
        multiAlignByLoc.append(multiAlignNoSep[seqName])
    multiAlignByLoc = zip(*multiAlignByLoc)
    return seqNameAll, multiAlignByLoc


def align_reduce(alignInSeg):
    """
    form alignment results of several segments to a single alignment
    only alignment of substitutions is used
    """
    res = []
    for aligni in alignInSeg.values():
        res += aligni
    return res


def align_str(alignReduced):
    """
    transform all alignments into two strings
    alignments involving '-' are removed
    """
    alignTransOnly = [ali for ali in alignReduced if ali.count('-') == 0]
    str1, str2 = zip(*alignTransOnly)
    str1 = ''.join(str1)
    str2 = ''.join(str2)
    return str1, str2


def align_pool(alignRes):
    vas = alignRes.values()
    res = []
    for va in vas:
        res = res + va
    return res


def lan_pair_list(lists):
    """
    construct a list of all language pairs, sorted
    """
    lanNames = lists.keys()
    lanNames.sort()
    lanPairsIter = combinations(lanNames, 2)
    lanPairsList = [lanPair for lanPair in lanPairsIter]
    return lanPairsList
