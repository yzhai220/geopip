# CTMC+NJ method.

import numpy as np
import scipy as sp

from io_json import *
from tree_util import get_out_name_dist_tree_files, get_rscript
from align_util import pair_align_from_multi_align_subs_only
from pip_est import opt_qmat_em_full, tree_use_r_for_unknown_number_of_leaves


# This function is not optimized for larger cList.
def nllk_2str_ctmc_bonly(b, psa, piProb, qMat, cList, qRates=[1]):
    """
    nllk of 2 strings based on CTMC model, for optimizing b only
    input:
        psa: pairwise string alignment, for example, [('A', 'T'), ('A', 'A'),
            ...]
        cList: a list of characters
        qRates: a list of qRates
    return:
        nllk
    """
    if b <= 0:
        # Set lower bound, and do not accept negative values.
        return 1.e10
    if b > 100:
        # Set upper bound, do accept bigger values, but truncate.
        return 1.e-10
    nChar = len(cList)
    pairList = [(ch1, ch2) for ch1 in cList for ch2 in cList]
    pairCount = [psa.count(pair) for pair in pairList]
    # startProb = np.repeat(piProb, nChar)
    tranProbAll = np.zeros(nChar * nChar)
    for qRate in qRates:
        tranProb = sp.linalg.expm(qMat * b * qRate)
        # print tranProb
        tranProb = tranProb.reshape(nChar * nChar)
        tranProbAll += tranProb
    tranProbAll = tranProbAll / len(qRates)
    # nllk = -np.sum((np.log(startProb) + np.log(tranProb)) * pairCount)
    # not that startProb is not influenced by b
    nllk = -np.sum(np.log(tranProbAll) * pairCount)
    return nllk


def opt_2str_ctmc_bonly(psa, piProb, qMat, cList, qRates=[1]):
    """optimization: for b only, with fixed rate matrix"""
    res = sp.optimize.minimize_scalar(nllk_2str_ctmc_bonly, args=(psa, piProb, qMat, cList, qRates))
    est = res.x
    nllk = res.fun
    return est, nllk


def opt_nstr_ctmc_bonly(pairAlignSubsOnly, piProb, qMat, cList, qRates=[1]):
    """
    optimization for all pairs, b only
    """
    m = {}    # record estimate
    nllkAll = 0    # record negative log-likelihood
    for pair, pairAlignAllSeg in pairAlignSubsOnly.iteritems():
        psa = pairAlignAllSeg[0]
        est, nllk = opt_2str_ctmc_bonly(psa, piProb, qMat, cList, qRates)
        m[pair] = est
        nllkAll += nllk
        print pair, '=', est
    return m, nllkAll


def opt_ctmc_full(qMat, multiAlign, javaDirectory, modelDirectory, eStepFile, parametersPath, inputLoc, outputLoc, dataLoc, execsLoc, rFileLoc, cList, qRates=[1], suffix='', updateQ=True, tol=1.e-2, bTol=1.e-3, iterMax=100):
    """
    optimization for all parameters in the CTMC model: qMat, tree (bDict)
    update qMat and tree (bDict) iteratively
    """
    outNameLoc, outDistLoc, outTreeLoc = get_out_name_dist_tree_files(dataLoc, suffix)
    rCodeNj = get_rscript(outNameLoc, outDistLoc, outTreeLoc, rFileLoc)
    print 'CTMC estimate for data in: %s' % (inputLoc)
    pairAlignSubsOnly = pair_align_from_multi_align_subs_only(multiAlign)
    piProb = pi_from_qmat(qMat)
    pairsList = pairAlignSubsOnly.keys()
    print '###  initialize tree  ###'
    bDict, nllkAll = opt_nstr_ctmc_bonly(pairAlignSubsOnly, piProb, qMat, cList, qRates)
    tree = tree_use_r_for_unknown_number_of_leaves(bDict, pairsList, rCodeNj, dataLoc, outTreeLoc, suffix, rooted=True)
    outAlignFile = inputLoc + '/' + 'all.align.txt'
    outTreeFile = inputLoc + '/' + 'all.tree.txt'
    dict_write_align_fasta(multiAlign, outAlignFile)
    write_tree(tree, outTreeFile)
    dif = 1.e10
    bDictRelativeDif = 1
    iterNum = 1
    while ((dif > tol) and (bDictRelativeDif > bTol) and (iterNum < iterMax)):
        if updateQ:
            print '### updating rate matrix Q ###\n'
            # write_align(pairsList, pairAlignSubsOnly, bDict, inputLoc)
            qMatNew, piProbNew = opt_qmat_em_full(qMat, cList, inputLoc, outputLoc, javaDirectory, modelDirectory, eStepFile, parametersPath, execsLoc)
            qMatRelativeDifMat = abs(qMatNew - qMat) / abs(qMat)
            qMatRelativeDif = qMatRelativeDifMat.max()
            qMat = qMatNew
            piProb = piProbNew
        else:
            print '### fixing rate matrix Q ###\n'
            qMatRelativeDif = 0
        print '### updating tree ###'
        bDictNew, nllkAll = opt_nstr_ctmc_bonly(pairAlignSubsOnly, piProb, qMat, cList, qRates)
        bDictRelativeDifVec = [abs(bDictNew[key] - bDict[key]) / bDict[key] for key in bDict.keys()]
        bDictRelativeDif = np.array(bDictRelativeDifVec).max()
        bDict = bDictNew
        print 'iter=%s: Q diff = %s, bDict diff = %s' % (iterNum, qMatRelativeDif, bDictRelativeDif)
        dif = max(qMatRelativeDif, bDictRelativeDif)
        tree = tree_use_r_for_unknown_number_of_leaves(bDict, pairsList, rCodeNj, dataLoc, outTreeLoc, suffix='', rooted=True)
        write_tree(tree, outTreeFile)
        iterNum += 1
    if iterNum == iterMax:
        print 'maximum iteration %d reached' % (iterMax)
    else:
        print 'optimization sucess!'
        print 'nllk = %s' % (nllkAll)
    tree = tree_use_r_for_unknown_number_of_leaves(bDict, pairsList, rCodeNj, dataLoc, outTreeLoc, suffix, rooted=True)
    write_tree(tree, outTreeFile)
    return qMat, bDict, tree, nllkAll
