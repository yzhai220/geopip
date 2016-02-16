# Inference of segmentation using dynamic programming for multiple alignments.

import numpy as np

from pip_msa import *


def logprob_align_cons_mat(msaList, ratesList, seqNames, tree, qMat, piProb, cList):
    """
    calculate the alignment log probability of for all single alignment under all possible rates
    calculate the constant part in the PIP density, in log scale
    input:
        msaList: a list of n multiple alignments, each of 'ATC' or ('A', 'T', 'C')
        ratesList: list of m possible rates
        seqNames: names of sequences
        tree:
        qMat, piProb: rate matrix and stationary distribution
        cList: list of basic characters, by default, [A, C, G, T] for DNA
    output:
        logProbAlignMat: matrix, n*m dim, log-probability of msa at one locus
        logProbConsMat: matrix, n*m dim, log-probability of psi(cPhi, m)
    """
    nRates = len(ratesList)
    nAligns = len(msaList)
    logProbAlignMat = np.empty((nAligns, nRates))
    logProbConsMat = np.empty((nAligns, nRates))
    nLeaf = len(seqNames)
    tau = tree.length()
    cPhi = '-' * nLeaf
    piProbExt = np.append(piProb, 0)
    cListExt = cList + ['-']
    msaUnique = set(msaList)
    for index in xrange(nRates):
        iRate, dRate = ratesList[index]
        # pre-calculate logProbAlignMat
        qMatExt = q_to_qext(qMat, dRate)
        logProbMsaUnique = logprob_msa_multi_site(msaUnique, seqNames, tree, qMatExt, piProbExt, dRate, cListExt)
        logProbMsaDict = dict(zip(msaUnique, logProbMsaUnique))
        # pre-calculate logProbConsMat
        nu = iRate * (tau + 1. / dRate)
        pc0 = prob_msa_one_site(cPhi, seqNames, tree, qMatExt, piProbExt, dRate, cListExt)
        logNu = np.log(nu)
        logProbPhi = (pc0 - 1) * nu
        for n in xrange(nAligns):
            msa = msaList[n]
            logProbAlignMat[n, index] = logProbMsaDict[msa]
            logProbPhi = logProbPhi - np.log(n+1) + logNu
            logProbConsMat[n, index] = logProbPhi
    return logProbAlignMat, logProbConsMat


def llh_max_geo_dp(alignLogProbMat, consLogProbMat, piProbRates, p, segMin=5):
    """
    calculate likelihood of BEST segmentation with geometric number of segments
    input:
        alignLogProbMat: a matrix of single alignment probabilities
        consLogProbMat: a matrix of constant part probability in PIP, related to alignment lengths as well
        piProbRates: stationary distributi of rate pairs
        p: the parameter of geometric distribution for segments
    output:
    """
    nAligns, nRates = consLogProbMat.shape
    f = np.zeros(nAligns)   # track segmentation of maximum likleihood
    g = np.zeros(nAligns)   # track the rate of last segment
    # log m_{i,j} + log pi_{j}
    consPiLogProbMat = consLogProbMat + np.log(piProbRates)
    llh = np.zeros(1)
    q = 1. - p
    logp = np.log(p)
    logq = np.log(q)
    for i in xrange(nAligns):
        llhNewAll = np.zeros((i+1, nRates))
        pLogVec = logp - logq
        for j in xrange(i+1):
            pLogVecNew = alignLogProbMat[i-j, :]
            pLogVec += pLogVecNew
            llhNewAll[i-j,:] = llh[i-j] + pLogVec + consPiLogProbMat[j, :] + logq
        llhNewAll[max(1, (i+2-segMin)):(i+1), :] = -1.e10
        llhNew = llhNewAll.max()
        index = llhNewAll.argmax()
        llh = np.append(llh, llhNew)
        f[i] = index / nRates
        g[i] = index % nRates
        f = f.astype(int)
        g = g.astype(int)
    return f, g


def get_seg_len_rate_from_f_g(f, g):
    """
    reconstruct the segmentation lengths and rates from the backward track function f and g
    input:
        f: from llh_max_geo function, a backward track function {0,...,n-1}->{0,...,n-1} tracking segment length
        g: from llh_max_geo function, a backward track function {0,...,n-1}->{0,...,m-1} tracking rates
    output:
        rateSegs: rates indices, a vectoer
    """
    fLen = len(f)
    fIndex = fLen - 1
    lenSegs = []
    rateSegs = []
    while fIndex>0:
        fIndexNew = f[fIndex] - 1
        lenSegNew = fIndex - fIndexNew
        rateSegNew = g[fIndex]
        lenSegs.append(lenSegNew)
        rateSegs.append(rateSegNew)
        fIndex = fIndexNew
    lenSegs.reverse()
    rateSegs.reverse()
    return lenSegs, rateSegs


def mle_seg_len_rate(p, msaList, ratesList, seqNames, tree, qMat, piProb, piProbRates, cList):
    """
    get the maximum likelihood segmentation using dynamic programming
    input:
        msaList: list, of MSAs
        piProbRates: probability of each rate pair
    output:
        lenSegs: list, number of alignments in each segment
        rateSegs: list, the rate index of each segment
    """
    logProbAlignMat, logProbConsMat = logprob_align_cons_mat(msaList, ratesList, seqNames, tree, qMat, piProb, cList)
    f, g = llh_max_geo_dp(logProbAlignMat, logProbConsMat, piProbRates, p)
    lenSegs, rateSegs = get_seg_len_rate_from_f_g(f, g)
    return lenSegs, rateSegs


def segratedict_from_ratesegs(rateSegs, ratesList):
    """
    construct segRateDict from rateSegs and ratesList
    input:
        rateSegs: list, index of rate of ratesList in each segs
        ratesList: all rate pairs
    output:
        segRateDict: dict, segId -> rate pair in this seg
    """
    segRateDict = {}
    nSegs = len(rateSegs)
    for segId in xrange(nSegs):
        segRateDict[segId] = ratesList[rateSegs[segId]]
    return segRateDict
