# GeoPIP+NJ method for DNA data.

import numpy as np
from scipy.optimize import minimize_scalar

from pip_msa import *
from geopip_msa import *
from pip_est import *

from tree_util import get_out_name_dist_tree_files, get_rscript
from pip_util import *
from align_util import multi_align_sep_by_seg, pair_align_from_multi_align, get_multi_align_in_all_seg


def update_segmentation_in_multi_align(multiAlign, lenSegs):
    multiAlignInSeg, segIds = multi_align_sep_by_seg(multiAlign, lenSegs)
    return multiAlignInSeg, segIds


def update_segmentation_in_align_in_seg(lenSegs, rateSegs, ratesList, multiAlign):
    """
    update pairwise alignments in segments under new segmentation
    note that because the number of segments changes, the segRateDict also changes, which also need to be updated
    input:
        lenSegs: length of segments
        rateSegs: estimated rate index in each segment
        ratesList: list, all rates
        multiAlign: multiple alignment, seqName -> string
    output:
        alignInSeg: dict, pair->segId->pairwise alignment (non empty pairs only, i.e., no ['-', '-'])
        segRateDict: dict, segId->segRate
    """
    alignInSeg = pair_align_from_multi_align(multiAlign, lenSegs)
    segRateDict = segratedict_from_ratesegs(rateSegs, ratesList)
    return alignInSeg, segRateDict


def update_segmentation_in_multi_aling_all_seg(multiAlign, lenSegs):
    """
    update multiple alignments in segments under new segmentation
    note that because the number of segments changes, the segRateDict also changes, which also need to be updated
    input:
        lenSegs: length of segments
        rateSegs: estimated rate index in each segment
        ratesList: list, all rates
        multiAlign: multiple alignment, seqName -> string
    output:
        multiAlignAllSeg: dict, pair->segId->pairwise alignment (non empty pairs only, i.e., no ['-', '-'])
    """
    multiAlignInSeg, segIds = multi_align_sep_by_seg(multiAlign, lenSegs)
    multiAlignAllSeg = get_multi_align_in_all_seg(multiAlignInSeg, segIds)
    return multiAlignAllSeg


def update_piprobrates(rateSegs, m, addOneSmoothing=True):
    """
    update piProbRates, the stationary distrubution of rate pairs
    naive Bayes estimate, counting how many appears in the segmentation result
    input:
        rateSegs: a list of 0, ..., m-1 indicate segment rates, from get_seg_len_rate_from_f_g
        m: number of all possible segments
        addOneSmoothing: bool, indicate whether add 1 smoothing is used here, by default is True
    output:
        piProbRates: probability for each segment rates
    """
    count = np.array([rateSegs.count(index) for index in xrange(m)])
    if addOneSmoothing:
        count = count + 1.
    piProbRates = count / count.sum()
    return piProbRates


def update_p(nSeg):
    """
    update p, the geometric parameter using p = 1 / E(nSeg)
    input:
        nSeg: number of segments
    output:
        p: the gememetric parameter
    """
    if nSeg == 1:
        p = 1. / nSeg - 0.1
    else:
        p = 1. / nSeg
    return p


def nllk_nlists_rate_cluster_inseg_only_msa(rate, segIdsInOneCluster, multiAlignAllSeg, tree, piProb, qMat, cList):
    """
    calculate nllk of n taxa, for optimization of (iRate, dRate) in segments of same rates
    note that only parts related to iRate, dRate is calculated
    input:
        segIdsInOneCluster: a list of segment IDs of this rate (iRate, dRate)
        multiAlignAllSeg: dict, segId -> seqName -> string in seq with segId
    ouput:
        nllk
    """
    nllk = 0
    if min(rate) <= 0:
        return 1.e10
    iRate = rate[0]
    dRate = rate[1]
    for segId in segIdsInOneCluster:
        multiAlign = multiAlignAllSeg[segId]
        # TODO: THIS PART CAN BE IMPROVED, BY CALCULATING dict outside the loop
        nllkNew = -logprob_msa(multiAlign, tree, qMat, piProb, iRate, dRate, cList)
        nllk += nllkNew
    # print param
    return nllk


def nllk_msa_geopip_final(ratesList, segRateDict, multiAlignAllSeg, tree, piProb, qMat, cList, piProbRates, p, rateSegs):
    """
    calculate the log-probability of MSA in all segments based on PIP
    """
    nllk = 0
    segIds = segRateDict.keys()
    rateAll = segRateDict.values()
    rateUniqueSet = set(rateAll)
    rateUnique = list(rateUniqueSet)
    for rate in rateUnique:
        segIdsInOneCluster = [segId for segId in segIds if segRateDict[segId] == rate]
        nllk += nllk_nlists_rate_cluster_inseg_only_msa(rate, segIdsInOneCluster, multiAlignAllSeg, tree, piProb, qMat, cList)
    # IMPROVE!!!
    q = 1 - p
    nllk += - (len(segRateDict) - 1) * np.log(p) - np.log(q)
    rateSegsCount = np.array([rateSegs.count(i) for i in xrange(len(ratesList))])
    nllk += - np.sum(np.log(piProbRates) * rateSegsCount)
    return nllk


# TODO: improve the efficiency of this function later.
def nllk_nlists_rate_cluster_inseg_only_msa_drate(dRate, segIdsInOneCluster, multiAlignAllSeg, tree, piProb, qMat, cList):
    """
    calculate nllk of n taxa, for optimization of (iRate, dRate) in segments of same rates
    iRate is calculated using MLE based on dRate
    input:
        segIdsInOneCluster: a list of segment IDs of this rate (iRate, dRate)
        multiAlignAllSeg: dict, segId -> seqName -> string in seq with segId
    ouput:
        nllk
    """
    nllk = 0
    if dRate <= 0:
        return 1.e10
    iRateVec = []
    for segId in segIdsInOneCluster:
        multiAlign = multiAlignAllSeg[segId]
        seqNames = multiAlign.keys()
        mLen = len(multiAlign.values()[0])
        pc0 = pc0_from_dRate_and_tree(dRate, seqNames, tree, qMat, piProb, cList)
        iRateNew = mle_irate_given_drate(dRate, tree, pc0, mLen)
        iRateVec.append(iRateNew)
    iRate = sum(iRateVec) / len(iRateVec)
    for segId in segIdsInOneCluster:
        multiAlign = multiAlignAllSeg[segId]
        # TODO: THIS PART CAN BE IMPROVED, BY CALCULATING dict outside the loop
        nllkNew = -logprob_msa(multiAlign, tree, qMat, piProb, iRate, dRate, cList)
        nllk += nllkNew
    # print param
    return nllk


def opt_nlists_rate_cluster_inseg_only_msa_drate(segIdsInOneCluster, multiAlignAllSeg, tree, piProb, qMat, cList):
    """
    optimization for estimating (iRate, dRate) in segments of same rates
    input:
        bDict: a dictionary of pairwise distances, pair -> distance
    output:
        est = (iRateEst, dRateEst)
    """
    res = minimize_scalar(nllk_nlists_rate_cluster_inseg_only_msa_drate, args=(segIdsInOneCluster, multiAlignAllSeg, tree, piProb, qMat, cList), tol=1.e-4)
    dRate = res.x
    if dRate <= 0:
        dRate = 0
    iRateVec = []
    for segId in segIdsInOneCluster:
        multiAlign = multiAlignAllSeg[segId]
        seqNames = multiAlign.keys()
        mLen = len(multiAlign.values()[0])
        pc0 = pc0_from_dRate_and_tree(dRate, seqNames, tree, qMat, piProb, cList)
        iRateNew = mle_irate_given_drate(dRate, tree, pc0, mLen)
        iRateVec.append(iRateNew)
    iRate = sum(iRateVec) / len(iRateVec)
    est = (iRate, dRate)
    return est


def opt_nlists_rate_all_cluster_inseg_only_msa_drate(ratesList, segRateDict, multiAlignAllSeg, tree, piProb, qMat, cList):
    """
    optimization for estimating (iRate, dRate) in all segments
    segments of the same rates are estimated together
    input:
        segRateDict: a dictionary of rates to start optimization, segId -> (iRateStart, dRateStart)
    output:
        res: a dictionary of est. rates in all segments, segId -> (iRateEst,
            dRateEst)
    """
    m = len(ratesList)
    est = {}
    segRateDictNew = {}
    segIds = segRateDict.keys()
    rateAll = segRateDict.values()
    rateUniqueSet = set(rateAll)
    rateUnique = list(rateUniqueSet)
    rateNotUsed = [rate for rate in ratesList if rate not in rateUnique]
    for rate in rateUnique:
        segIdsInOneCluster = [segId for segId in segIds if segRateDict[segId] == rate]
        tem = opt_nlists_rate_cluster_inseg_only_msa_drate(segIdsInOneCluster, multiAlignAllSeg, tree, piProb, qMat, cList)
        est[rate] = tem
    for segId in segIds:
        segRateDictNew[segId] = est[segRateDict[segId]]
    # update ratesList
    ratesListNew = []
    for index in xrange(m):
        rate = ratesList[index]
        if rate in rateUnique:
            ratesListNew += [est[rate]]
        elif rate in rateNotUsed:
            ratesListNew += [rate]
        else:
            print 'Error: rate appears but not updated!'
    ratesListNew.sort()
    return segRateDictNew, ratesListNew


# TODO: improve the efficiency of this function later.
def nllk_nlists_rate_all_cluster_inseg_only_msa_drate_fix_drate_times_b(delta, ratesList, segRateDict, multiAlignAllSeg, tree, piProb, qMat, cList):
    """
    nllk of all segments
    input:
        segRateDict: a dictionary of rates to start optimization, segId -> (iRateStart, dRateStart)
    output:
        res: a dictionary of est. rates in all segments, segId -> (iRateEst,
            dRateEst)
    """
    if delta <= 0:
        return 1.e10
    nllk = 0
    segIds = segRateDict.keys()
    rateAll = segRateDict.values()
    rateUniqueSet = set(rateAll)
    rateUnique = list(rateUniqueSet)
    treeScaled = deepcopy(tree)
    treeScaled.scale_edges(1./delta)
    for rate in rateUnique:
        dRate = rate[1] * delta
        segIdsInOneCluster = [segId for segId in segIds if segRateDict[segId] == rate]
        nllk += nllk_nlists_rate_cluster_inseg_only_msa_drate(dRate, segIdsInOneCluster, multiAlignAllSeg, treeScaled, piProb, qMat, cList)
    # update ratesList
    return nllk


def opt_nlists_rate_all_cluster_inseg_only_msa_drate_fix_drate_times_b(ratesList, segRateDict, multiAlignAllSeg, tree, piProb, qMat, cList):
    """
    optimization for estimating (iRate, dRate) in all segments when dRate*tau is fixed
    """
    res = minimize_scalar(nllk_nlists_rate_all_cluster_inseg_only_msa_drate_fix_drate_times_b, args=(ratesList, segRateDict, multiAlignAllSeg, tree, piProb, qMat, cList), tol=1.e-4)
    delta = res.x
    print delta
    treeNew = deepcopy(tree)
    treeNew.scale_edges(1/delta)
    segIds = segRateDict.keys()
    rateAll = segRateDict.values()
    rateUniqueSet = set(rateAll)
    rateUnique = list(rateUniqueSet)
    segRateDictNew = deepcopy(segRateDict)
    ratesListNew = deepcopy(ratesList)
    for rate in rateUnique:
        iRateVec = []
        dRate = rate[1] * delta
        segIdsInOneCluster = [segId for segId in segIds if segRateDict[segId] == rate]
        for segId in segIdsInOneCluster:
            multiAlign = multiAlignAllSeg[segId]
            seqNames = multiAlign.keys()
            mLen = len(multiAlign.values()[0])
            pc0 = pc0_from_dRate_and_tree(dRate, seqNames, treeNew, qMat, piProb, cList)
            iRateNew = mle_irate_given_drate(dRate, treeNew, pc0, mLen)
            iRateVec.append(iRateNew)
        iRate = sum(iRateVec) / len(iRateVec)
        rateNew = (iRate, dRate)
        ratesListNew[ratesListNew.index(rate)] = rateNew
        for k, v in segRateDictNew.iteritems():
            if v == rate:
                segRateDictNew[k] = rateNew
            else:
                segRateDictNew[k] = v
    return ratesListNew, segRateDictNew, treeNew


def opt_geopip_full(m, p, qMat, segRateDict, piProbRates, ratesList, multiAlign, lenSegs, javaDirectory, modelDirectory, eStepFile, parametersPath, inputLoc, outputLoc, dataLoc, execsLoc, rFileLoc, cList, qRates=[1.], suffix='', updateQ=True, updateSeg=True, updateRate=True, updateRateFixdRateTimesTau=True, rooted=True, tol=1.e-2, bTol=1.e-3, iterMax=100):
    """
    main function for optimization
    updating in segment rates, out segment rates, qMat, bDict iteratively
    we estimate bDict first, so bDict is not needed as input
    """
    outNameLoc, outDistLoc, outTreeLoc = get_out_name_dist_tree_files(dataLoc, suffix)
    rCodeNj = get_rscript(outNameLoc, outDistLoc, outTreeLoc, rFileLoc)
    seqNames = multiAlign.keys()
    msaList = zip(*multiAlign.values())
    print 'simulation run: %s' % (inputLoc)
    multiAlignAllSeg = update_segmentation_in_multi_aling_all_seg(multiAlign, lenSegs)
    alignsInSeg = pair_align_from_multi_align(multiAlign, lenSegs)
    pairsList = alignsInSeg.keys()
    piProb = pi_from_qmat(qMat)
    print '###  initialize branch lengths  ###'
    bDict = opt_nlists_bonly(pairsList, alignsInSeg, segRateDict, piProb, qMat, cList, qRates)
    # write_dist_from_bdict(bDict, dataLoc, suffix)    # this function can be improved, by using outNameLoc, outDistLoc, outTreeLoc instead
    tree = tree_use_r_for_unknown_number_of_leaves(bDict, pairsList, rCodeNj, dataLoc, outTreeLoc, suffix, rooted)
    outAlignFile = inputLoc + '/' + 'all.align.txt'
    outTreeFile = inputLoc + '/' + 'all.tree.txt'
    dict_write_align_fasta(multiAlign, outAlignFile)
    write_tree(tree, outTreeFile)
    tree.reroot_at_midpoint()
    dif = 1.e10
    bDictRelativeDif = 1
    iterNum = 1
    nllk = 1.e10
    # nllkDif = 1
    while ((dif > tol) and (bDictRelativeDif > bTol) and (iterNum < iterMax)):
        if updateSeg:
            print '### update segmentation ###'
            lenSegs, rateSegs = mle_seg_len_rate(p, msaList, ratesList, seqNames, tree, qMat, piProb, piProbRates, cList)
            nSeg = len(rateSegs)
            alignsInSeg, segRateDict = update_segmentation_in_align_in_seg(lenSegs, rateSegs, ratesList, multiAlign)
            multiAlignAllSeg = update_segmentation_in_multi_aling_all_seg(multiAlign, lenSegs)
            piProbRates = update_piprobrates(rateSegs, m, True)
            p = update_p(nSeg)
        print lenSegs
        if updateRateFixdRateTimesTau:
            print '### updating insertion and deletion rate when dRate*tau is fixed ###'
            ratesListNew, segRateDictNew, treeNew = opt_nlists_rate_all_cluster_inseg_only_msa_drate_fix_drate_times_b(ratesList, segRateDict, multiAlignAllSeg, tree, piProb, qMat, cList)
            print 'rate:', ratesList
            print 'new rate:', ratesListNew
            dRateRelativeDifFixdRateTimesTauList = [abs(ratesListNew[index][1] - ratesList[index][1]) / ratesList[index][1] for index in xrange(m) if ratesListNew[index][1] > 1.e-3]
            dRateRelativeDifFixdRateTimesTau = np.array(dRateRelativeDifFixdRateTimesTauList).max()
            ratesList = ratesListNew
            segRateDict = segRateDictNew
            tree = treeNew
        else:
            dRateRelativeDifFixdRateTimesTau = 0
        if updateRate:
            print '### updating in segment rates ###'
            # segRateDictNew, ratesListNew = opt_nlists_rate_all_cluster_inseg_only_msa(ratesList, segRateDict, multiAlignAllSeg, tree, piProb, qMat, cList)
            segRateDictNew, ratesListNew = opt_nlists_rate_all_cluster_inseg_only_msa_drate(ratesList, segRateDict, multiAlignAllSeg, tree, piProb, qMat, cList)
            print 'rate:', ratesList
            print 'new rate:', ratesListNew
            dRateRelativeDifList = [abs(ratesListNew[index][1] - ratesList[index][1]) / ratesList[index][1] for index in xrange(m) if ratesListNew[index][1] > 1.e-3]
            dRateRelativeDif = np.array(dRateRelativeDifList).max()
            segRateDict = segRateDictNew
            ratesList = ratesListNew
        else:
            dRateRelativeDif = 0
            print 'rate is Fixed:', ratesList
        #####################
        # another approach based on clustering
        # segRateDictNew = opt_nlists_rate_all_inseg_only(segRateDict, alignsInSeg, bDict, pairsList, piProb, qMat, cList)
        # # IMPROVE THIS LATER, scale each rate separately first maybe?
        # inSegRateDif = max([abs(segRateDictNew[segId] - np.array(segRateDict[segId])).max() for segId in segRateDictNew.keys()])
        # # update ratesList
        # ratesList, segRateDictNew = ratesList_kmean(segRateDictNew, m)
        # print segRateDict
        # print segRateDictNew
        # segRateDict = segRateDictNew
        # display estiamteed in segment rates
        # for segId, segRate in segRateDict.iteritems():
        #     print 'seg id = %s, \t segRate = %s' %(segId, segRate)
        ######################
        if updateQ:
            print '### updating rate matrix Q ###'
            # if alignUpdate:
            #     aligns = align_mle(lans, pairsList, bDict, iRate, dRate, qMat, piProb, cList)
            # if cogUpdate:
            #     alignsCogOnlyEst = est_align_cog(lans, pairsList, iRate, dRate, cRate, cRateVec, bDict, piProb, qMat, cList, aligns, criticalValue)
            # else:
            #     alignsCogOnlyEst = alignsCogOnly
            # write_align(pairsList, alignsInSeg, bDict, inputLoc)
            qMatNew, piProbNew = opt_qmat_em_full(qMat, cList, inputLoc, outputLoc, javaDirectory, modelDirectory, eStepFile, parametersPath, execsLoc, lEstepTol=1.e-2, iterMax=100)
            qMatRelativeDifMat = abs(qMatNew - qMat) / abs(qMat)
            qMatRelativeDif = qMatRelativeDifMat.max()
            # qMatDif = abs(qMatNew - qMat).max()
            qMat = qMatNew
            piProb = piProbNew
            # qMatDif = 0
        else:
            print '### fixing rate matrix Q ###'
            qMatRelativeDif = 0
        print '### updating branch lengths ###'
        bDictNew = opt_nlists_bonly(pairsList, alignsInSeg, segRateDict, piProb, qMat, cList, qRates)
        # bDictDifVec = [abs(bDictNew[key] - bDict[key]) for key in bDict.keys()]
        bDictRelativeDifVec = [abs(bDictNew[key] - bDict[key]) / bDict[key] for key in bDict.keys()]
        bDictRelativeDif = np.array(bDictRelativeDifVec).max()
        # bDictDif = np.array(bDictDifVec).max()
        bDict = bDictNew
        # write_dist_from_bdict(bDict, dataLoc, suffix)
        tree = tree_use_r_for_unknown_number_of_leaves(bDict, pairsList, rCodeNj, dataLoc, outTreeLoc, suffix, rooted)
        print 'iter=%s: rates in seg diff = %s, rates in seg diff (fix dRate*tau) = %s, Q diff = %s, b diff = %s' % (iterNum, dRateRelativeDif, dRateRelativeDifFixdRateTimesTau, qMatRelativeDif, bDictRelativeDif)
        # print 'iter=%s: rates in seg diff = %s, Q diff = %s, b diff = %s' % (iterNum, segRateDif, qMatDif, bDictDif)
        dif = max(dRateRelativeDif, dRateRelativeDifFixdRateTimesTau, qMatRelativeDif, bDictRelativeDif)
        tree.reroot_at_midpoint()
        nllkNew = nllk_msa_geopip_final(ratesList, segRateDict, multiAlignAllSeg, tree, piProb, qMat, cList, piProbRates, p, rateSegs)
        nllkDif = -nllkNew + nllk
        nllkOld = nllk
        nllk = nllkNew
        print 'llk increase =', nllkDif
        if nllkDif <= 0:
            print 'Log-Likelihood is decreasing! BREAK!'
            bDict = bDictOld
            qMat = qMatOld
            segRateDict = segRateDictOld
            alignsInSeg = alignsInSegOld
            p = pOld
            piProbRates = piProbRatesOld
            ratesList = ratesListOld
            lenSegs = lenSegsOld
            tree = treeOld
            nllk = nllkOld
            break
        # do not update fixing dRate*tau if changes is small
        if dRateRelativeDifFixdRateTimesTau < 0.05:
            updateRateFixdRateTimesTau = False
        write_tree(tree, outTreeFile)
        iterNum += 1
        bDictOld = bDict
        qMatOld = qMat
        segRateDictOld = segRateDict
        alignsInSegOld = alignsInSeg
        pOld = p
        piProbRatesOld = piProbRates
        ratesListOld = ratesList
        lenSegsOld = lenSegs
        treeOld = tree
    if iterNum == iterMax:
        print 'maximum iteration %d reached' % (iterMax)
    else:
        print 'optimization sucess!'
    return bDict, qMat, segRateDict, alignsInSeg, p, piProbRates, ratesList, lenSegs, tree, nllk
