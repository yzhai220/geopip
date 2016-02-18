# Estimate parameters in PIP.
# Parameters to be estimated: iRate, dRate, qMat, tree (bDict).

from scipy.optimize import minimize, minimize_scalar
import dendropy
import copy

from io_json import *
from align_util import *
from pip_util import pi_from_qmat, q_to_qext
from pip_msa import logprob_align, logprob_msa, prob_msa_one_site
from tree_util import get_out_name_dist_tree_files, get_rscript, nj_tree_from_bdict_using_r


def mle_irate_given_drate(dRate, tree, pc0, mLen):
    """
    For any given dRate, then iRate can be calculated using MLE
    """
    tau = tree.length()
    nu = mLen / (1. - pc0)
    iRate = nu / (tau + 1. / dRate)
    return iRate


def pc0_from_dRate_and_tree(dRate, seqNames, tree, qMat, piProb, cList, qRates=[1.]):
    """
    calculate P(cPhi) under PIP
    """
    nLeaf = len(seqNames)
    cPhi = '-' * nLeaf
    cListExt = cList + ['-']
    pc0 = 0
    for qRate in qRates:
        qMatExt = q_to_qext(qMat*qRate, dRate)
        piProbExt = np.append(piProb, 0)
        pc0 += prob_msa_one_site(cPhi, seqNames, tree, qMatExt, piProbExt, dRate, cListExt)
    pc0 = pc0 / len(qRates)
    return pc0


def align_inseg_logprob(alignRes, segRateDict, piProb, qMat, b, cList, qRates=[1.]):
    """
    calculate the probability of alignments of sequences within each segment
    input:
        alignRes: alignemnt of sequences wthin EACH segments
        segRateDict: dictionary of segment Id to rate=(iRate, dRate) in this segment
        piProb: stationary distribution of qMat
        qMat: rate matrix
        b: branch length
        cList: list of basic characters, by default, [A, C, G, T] for DNA
        qRates: list of rates facters for substitution rate variation.
    output:
        probAlign: a dictionary of segment ID to likelihood
    """
    logProbAlign = {}
    for segId, segValue in alignRes.iteritems():
        rate = segRateDict[segId]
        iRate = rate[0]
        dRate = rate[1]
        logProbAlign[segId] = logprob_align(b, segValue, piProb, qMat, iRate, dRate, cList, qRates)
    return logProbAlign


def nllk_2lists_bonly(b, alignInSeg, segRateDict, piProb, qMat, cList, qRates=[1.]):
    """
    geo PIP model
    nllk for 2 lists: used to estimate b only
    input:
        alignInSeq: alignment of sequences within each segment, a dictionary
        ratesList: character list of segments, specified by (iRate, dRate)
        segRateDic: dictionary of segment id to within segmant (iRate, dRate)
        piProb, qMat, stationary distribution and rate matrix of basic elements
        cList: character list of basic elements, for DNA, it is [A, C, G, T]
        qRates: list of rates facters for substitution rate variation.
    output:
        nllk: negative log likelihood
    """
    if b <= 0:    # set lower bound, do not accept negative values
        return 1.e10
    if b > 100:    # set upper bound, do accept bigger values, but truncate
        return 1.e-10    # if bigger than 100, then return 100
    nllk = 0
    for qRate in qRates:
        bScaled = b * qRate
        logProbAlignInSeg = align_inseg_logprob(alignInSeg, segRateDict, piProb, qMat, bScaled, cList, qRates)
        nllk += - np.sum(logProbAlignInSeg.values())
    return nllk


def opt_2lists_bonly(alignInSeg, segRateDict, piProb, qMat, cList, qRates=[1.]):
    """
    geo PIP model
    optimization for branch length between 2 sequences
    input: same as in nllk_2lists_bonly
    output: MLE of branch length
    """
    res = minimize_scalar(nllk_2lists_bonly, args=(alignInSeg, segRateDict, piProb, qMat, cList, qRates))
    est = res.x
    return est


def opt_nlists_bonly(pairsList, alignsInSeg, segRateDict, piProb, qMat, cList, qRates=[1.]):
    """
    optimization for n lists, estimating pairwise branch lengths only
    input:
        pairsList: list of all pairs of taxa
        alignsInSeg: alignment of sequences within each segments, a dictionary, pair -> alignInSeg
    output:
        m: a dictionary of estimated pairwise distances, pair -> est. dist.
    """
    m = {}    # record estimate
    for pair in pairsList:
        alignInSeg = alignsInSeg[pair]
        est = opt_2lists_bonly(alignInSeg, segRateDict, piProb, qMat, cList, qRates)
        m[pair] = est
        print pair, '=', est
    return m


def tree_use_r_for_unknown_number_of_leaves(bDict, pairsList, rCodeNj, dataLoc, outTreeLoc, suffix='', rooted=True):
    if len(bDict) == 1:
    # only two leaves, do not run R code, just mid-point rooting
        lan1 = pairsList[0][0]
        lan2 = pairsList[0][1]
        b = bDict.values()[0]
        treeStr = '(' + str(lan1) + ':' + str(b/2) + ',' + str(lan2) + ':' + str(b/2) + ');'
        tree = dendropy.Tree.get_from_string(treeStr, schema="newick")
    else:
    # three or more leaves, run R code of NJ and mid-point rooting
        tree = nj_tree_from_bdict_using_r(rCodeNj, bDict, dataLoc, outTreeLoc, suffix, rooted=True)
    return tree


def nllk_rate(rate, multiAlign, tree, qMat, piProb, cList):
    """
    negative log-likelihood for rate = (iRate, dRate) based on PIP
    """
    iRate, dRate = rate
    if iRate <= 0 or dRate <= 0:
        return 1.e10
    nllk = - logprob_msa(multiAlign, tree, qMat, piProb, iRate, dRate, cList)
    return nllk


def opt_rate(rateStart, multiAlign, tree, qMat, piProb, cList):
    """
    optimization function for rate = (iRate, dRate) based on PIP
    """
    res = minimize(nllk_rate, rateStart, method='BFGS', options={'disp': True}, tol=1.e-2, args=(multiAlign, tree, qMat, piProb, cList))
    est = tuple(res.x)
    return est


def nllk_drate(dRate, multiAlign, tree, qMat, piProb, cList, qRates=[1.]):
    """
    negative log-likelihood for rate = (iRate, dRate) based on PIP
    """
    if dRate <= 0:
        return 1.e10
    seqNames = multiAlign.keys()
    mLen = len(multiAlign.values()[0])
    pc0 = pc0_from_dRate_and_tree(dRate, seqNames, tree, qMat, piProb, cList, qRates)
    iRate = mle_irate_given_drate(dRate, tree, pc0, mLen)
    nllk = - logprob_msa(multiAlign, tree, qMat, piProb, iRate, dRate, cList, qRates)
    return nllk


def opt_drate(multiAlign, tree, qMat, piProb, cList, qRates=[1.]):
    """
    optimization function for rate = (iRate, dRate) based on PIP
    """
    res = minimize_scalar(nllk_drate, args=(multiAlign, tree, qMat, piProb, cList, qRates))
    dRate = res.x
    seqNames = multiAlign.keys()
    nLeaf = len(seqNames)
    cPhi = '-' * nLeaf
    tau = tree.length()
    cListExt = cList + ['-']
    piProbExt = np.append(piProb, 0)
    pc0 = 0
    for qRate in qRates:
        qMatExt = q_to_qext(qMat*qRate, dRate)
        pc0 += prob_msa_one_site(cPhi, seqNames, tree, qMatExt, piProbExt, dRate, cListExt)
    pc0 = pc0 / len(qRates)
    mlen = len(multiAlign.values()[0])
    nu = mlen / (1. - pc0)
    # logPsi = -np.sum(np.log(np.arange(1, mlen+1))) + mlen * np.log(nu) + (pc0 - 1) * nu
    # nu = iRate * (tau + 1. / dRate)
    iRate = nu / (tau + 1. / dRate)
    rate = (iRate, dRate)
    return rate


def opt_qmat_em_full(qMat, cList, inputLoc, outputLoc, javaDirectory, modelDirectory, eStepFile, parametersPath, execsLoc, lEstepTol=1.e-2, iterMax=100):
    """
    optimization function combining E step amd M step
    NEED GLOBAL VARIABLES: parametersPath, inputLoc, outputLoc, eStepFile, execsLoc
    lEstepTol: tolerance for likelihood in the Estep
    CURRENTLY: ONLY MONITOR LIKELIHOOD
    TO ADD LATER: qTol, wTol: tolerance for qMat, weight
    """
    nCharType = len(cList)
    treeAlignDict = get_treefiles_alignfiles(inputLoc)
    javacodeEBatch = java_estep_batch(treeAlignDict, outputLoc, parametersPath, javaDirectory)
    # execFile = find_newest_exec(execsLoc)
    javacodeM = java_mstep(modelDirectory, eStepFile, javaDirectory)
    # prepare for E step
    write_param_for_estep(qMat, parametersPath)
    # E step: first run
    [os.system(javacodeE) for javacodeE in javacodeEBatch]
    # collect all file names needed from first E step run
    listFiles, rootFiles, llhFiles = filepath_from_output_folder(outputLoc)
    llh = llh_from_llhfiles(llhFiles)
    lEstepDiff = 1.0
    iterCount = 0
    while (lEstepDiff > lEstepTol) and (iterCount < iterMax):
        # prepare for M step
        write_estepfiles_for_mstep(listFiles, rootFiles, eStepFile, nCharType, cList)
        # M Step
        os.system(javacodeM)
        execFile = find_newest_exec(execsLoc)
        qMatNew, piProbNew = write_param_estep_from_mstep(execFile, parametersPath)
        # E Step
        [os.system(javacodeE) for javacodeE in javacodeEBatch]
        llhNew = llh_from_llhfiles(llhFiles)
        lEstepDiff = abs(llhNew - llh)
        llh = llhNew
        iterCount += 1
        print '################################'
        print 'llhDiff = %s, \t EM iter = %s \n' % (lEstepDiff, iterCount)
    return qMatNew, piProbNew


def opt_pip_full(rate, qMat, multiAlign, javaDirectory, modelDirectory, eStepFile, parametersPath, inputLoc, outputLoc, dataLoc, execsLoc, rFileLoc, cList, qRates=[1.], suffix='', updateQ=True, updateRate=True, updateRateFixdRateTimesb=True, tol=1.e-2, bTol=1.e-3, iterMax=100):
    """
    optimization for all parameters: iRate, dRate, qMat, tree (bDict) in PIP
    updating in iRate and dRate, qMat, tree (bDict) iteratively
    estimate bDict first given other, so starting bDict (tree) is not needed
    """
    outNameLoc, outDistLoc, outTreeLoc = get_out_name_dist_tree_files(dataLoc, suffix)
    rCodeNj = get_rscript(outNameLoc, outDistLoc, outTreeLoc, rFileLoc)
    print 'simulation run: %s' % (inputLoc)
    # directory for runing EM
    alignInSeg = pair_align_from_multi_align(multiAlign)
    pairsList = alignInSeg.keys()
    piProb = pi_from_qmat(qMat)
    print '###  initialize tree  ###'
    segRateDict = {0: rate}   # only one segment
    bDict = opt_nlists_bonly(pairsList, alignInSeg, segRateDict, piProb, qMat, cList, qRates)
    # write_dist_from_bdict(bDict, dataLoc, suffix)
    tree = tree_use_r_for_unknown_number_of_leaves(bDict, pairsList, rCodeNj, dataLoc, outTreeLoc, suffix, rooted=True)
    outAlignFile = inputLoc + '/' + 'all.align.txt'
    outTreeFile = inputLoc + '/' + 'all.tree.txt'
    dict_write_align_fasta(multiAlign, outAlignFile)
    write_tree(tree, outTreeFile)
    tree.reroot_at_midpoint()
    dif = 1.e10
    bDictRelativeDif = 1
    iterNum = 1
    nllk = 1.e10
    nllkDif = 1
    while ((dif > tol) and (bDictRelativeDif > bTol) and (iterNum < iterMax)):
        if updateRateFixdRateTimesb:
            print '### updating insertion and deletion rate when dRate*b is fixed ###'
            dRate = rate[1]
            print tree.length()
            rateNew, tree = opt_drate_fix_drate_times_b(multiAlign, dRate, tree, qMat, piProb, cList, qRates=[1.])
            print tree.length()
            print rate
            print rateNew
            print (np.array(rateNew) - np.array(rate)) / np.array(rate)
            dRateRelativeDifFixdRateTimesb = abs(rate[1] - rateNew[1]) / rate[1]
            rate = rateNew
        else:
            dRateRelativeDifFixdRateTimesb = 0
        if updateRate:
            print '### updating insertion and deletion rate ###\n'
            # rateNew = opt_rate(rate, multiAlign, tree, qMat, piProb, cList)
            rateNew = opt_drate(multiAlign, tree, qMat, piProb, cList, qRates)
            print rate
            print rateNew
            print (np.array(rateNew) - np.array(rate)) / np.array(rate)
            dRateRelativeDif = abs(rate[1] - rateNew[1]) / rate[1]
            rate = rateNew
        else:
            print 'rate is fixed:', rate
            dRateRelativeDif = 0
        if updateQ:
            print '### updating rate matrix Q ###\n'
            qMatNew, piProbNew = opt_qmat_em_full(qMat, cList, inputLoc, outputLoc, javaDirectory, modelDirectory, eStepFile, parametersPath, execsLoc)
            qMatRelativeDifMat = abs(qMatNew - qMat) / abs(qMat)
            qMatRelativeDif = qMatRelativeDifMat.max()
            qMat = qMatNew
            piProb = piProbNew
        else:
            print '### fixing rate matrix Q ###\n'
            qMatRelativeDif = 0
        print '### updating tree ###'
        segRateDict = {0: rate}   # only one segment
        bDictNew = opt_nlists_bonly(pairsList, alignInSeg, segRateDict, piProb, qMat, cList, qRates)
        bDictRelativeDifVec = [abs(bDictNew[key] - bDict[key]) / bDict[key] for key in bDict.keys()]
        bDictRelativeDif = np.array(bDictRelativeDifVec).max()
        print (np.array(bDictNew.values()) - np.array(bDict.values())) / np.array(bDict.values())
        bDict = bDictNew
        # write_dist_from_bdict(bDict, dataLoc, suffix)
        tree = tree_use_r_for_unknown_number_of_leaves(bDict, pairsList, rCodeNj, dataLoc, outTreeLoc, suffix, rooted=True)
        # print 'iter=%s: iRate diff = %s, Q diff = %s, bDict diff = %s' % (iterNum, iRateRelativeDif, qMatRelativeDif, bDictRelativeDif)
        # dif = max(iRateRelativeDif, qMatRelativeDif, bDictRelativeDif)
        print 'iter=%s: dRate diff = %s, dRate diff (fix dRate*b) = %s, Q diff = %s, bDict diff = %s' % (iterNum, dRateRelativeDif, dRateRelativeDifFixdRateTimesb, qMatRelativeDif, bDictRelativeDif)
        dif = max(dRateRelativeDif, qMatRelativeDif, bDictRelativeDif, dRateRelativeDifFixdRateTimesb)
        tree.reroot_at_midpoint()
        nllkNew = nllk_rate(rate, multiAlign, tree, qMat, piProb, cList)
        nllkDif = nllk - nllkNew
        print 'llk increase =', nllkDif
        nllk = nllkNew
        write_tree(tree, outTreeFile)
        iterNum += 1
        # if dRateRelativeDifFixdRateTimesb is small, then skip that update
        if dRateRelativeDifFixdRateTimesb < 0.5:
            updateRateFixdRateTimesb = False
        if nllkDif <= 0:
            print 'Log-Lilkelihood is decreasing! BREAK'
            break
    if iterNum == iterMax:
        print 'maximum iteration %d reached' % (iterMax)
    else:
        print 'optimization sucess!'
    rate = [rate]
    return rate, qMat, bDict, tree


def pip_start_data(multiAlign):
    """
    start value of (iRate, dRate) for inference using PIP based on data
    """
    values = multiAlign.values()
    lenMSA = len(values[0])
    nonEmtpyValueVec = [lenMSA - value.count('-') for value in values]
    lenNonEmptyMean = np.array(nonEmtpyValueVec).mean()
    dRate = 1. - lenNonEmptyMean / lenMSA
    iRate = dRate * lenNonEmptyMean
    rate = (iRate, dRate)
    return rate


def nllk_drate_fix_drate_times_b(dRate, multiAlign, dRateTimesTreeTotalBranchLength, tree, qMat, piProb, cList, qRates=[1.]):
    """
    negative log-likelihood for rate = (iRate, dRate) based on PIP
    """
    if dRate <= 0:
        return 1.e10
    seqNames = multiAlign.keys()
    mLen = len(multiAlign.values()[0])
    treeScaled = copy.deepcopy(tree)
    treeScaled.scale_edges(dRateTimesTreeTotalBranchLength / (dRate * tree.length()))
    # print treeScaled.length()
    # print dRate
    # print treeScaled.length() * dRate
    pc0 = pc0_from_dRate_and_tree(dRate, seqNames, treeScaled, qMat, piProb, cList, qRates)
    iRate = mle_irate_given_drate(dRate, treeScaled, pc0, mLen)
    nllk = - logprob_msa(multiAlign, treeScaled, qMat, piProb, iRate, dRate, cList, qRates)
    # print nllk
    return nllk


def opt_drate_fix_drate_times_b(multiAlign, dRate, tree, qMat, piProb, cList, qRates=[1.]):
    """
    optimization function for rate = (iRate, dRate) based on PIP
    """
    dRateTimesTreeTotalBranchLength = dRate * tree.length()
    res = minimize_scalar(nllk_drate_fix_drate_times_b, args=(multiAlign, dRateTimesTreeTotalBranchLength, tree, qMat, piProb, cList, qRates))
    dRate = res.x
    treeNew = copy.deepcopy(tree)
    treeNew.scale_edges(dRateTimesTreeTotalBranchLength / (dRate * treeNew.length()))
    seqNames = multiAlign.keys()
    nLeaf = len(seqNames)
    cPhi = '-' * nLeaf
    tau = tree.length()
    cListExt = cList + ['-']
    piProbExt = np.append(piProb, 0)
    pc0 = 0
    for qRate in qRates:
        qMatExt = q_to_qext(qMat*qRate, dRate)
        pc0 += prob_msa_one_site(cPhi, seqNames, tree, qMatExt, piProbExt, dRate, cListExt)
    pc0 = pc0 / len(qRates)
    mlen = len(multiAlign.values()[0])
    nu = mlen / (1. - pc0)
    # logPsi = -np.sum(np.log(np.arange(1, mlen+1))) + mlen * np.log(nu) + (pc0 - 1) * nu
    # nu = iRate * (tau + 1. / dRate)
    iRate = nu / (tau + 1. / dRate)
    rate = (iRate, dRate)
    return rate, treeNew
