# PIP probability for multiple string alignments.
# See Bouchard and Jordan (2013).
# "Evolutionary inference via the Poisson Indel Process".
# All equation references point to this paper.

import numpy as np
import scipy as sp
from copy import deepcopy

from pip_util import q_to_qext


def fv(c1, c2, piExt, pMatExt, b, dRate, cListExt):
    """
    function used to calculate PIP density
    suppose that p0 = p1, so b1 = 0, b2 = b  !!!
    c1, c2: alignment element in p1 and p2
    piExt: extended stationary distribution with 0 for '-'
    qMatExt: extended rate matrix, with deletion row and column
    b: branch length from p1 (p0) to p2
    dRate: deletion rate
    cListExt: character list, extended with '-'
    """
    index1 = cListExt.index(c1)
    index2 = cListExt.index(c2)
    pic1 = piExt[index1]
    pic2 = piExt[index2]
    logic1 = (c1 != '-')
    logic2 = (c2 != '-')
    if logic1:
        ft0 = pic1 * pMatExt[index1, index2]
    else:
        ft0 = 0
    f0 = ft0
    beta = (1 - np.exp(-dRate * b)) / (b * dRate)
    if ((not logic1) and logic2):
        f2 = beta * pic2
    elif ((not logic1) and (not logic2)):
        f2 = 1 - beta
    else:
        f2 = 0
    if (logic1 or logic2):   # not really necesarry
        f1 = 0   # not really necesarry
    else:   # not really necesarry
        f1 = 1   # not really necesarry
    return f1, f2, f0


def pv_pri(b1, b2, dRate):
    """
    prior of P(V = v)
    """
    tau = b1 + b2
    denom = tau + 1. / dRate
    p1 = b1 / denom
    p2 = b2 / denom
    p0 = 1. / dRate / denom
    return p1, p2, p0


def prob_c_all(piExt, pMat, b, dRate, cListExt):
    """
    Pr(C=c), i.e. p(c), the probability of one alignment for one location
    calculating all possible c's and store it for later
    """
    pc = {}
    pvPrior = pv_pri(0, b, dRate)
    for c1 in cListExt:
        pc[c1] = {}
        for c2 in cListExt:
            fvProb = fv(c1, c2, piExt, pMat, b, dRate, cListExt)
            pc[c1][c2] = np.dot(pvPrior, fvProb)
    return pc


def logprob_m(m, pc0, pcAll, b, iRate, dRate):
    """
    calculating log(p(m)), the log probability of a fixed alignment for more than one location
    input:
        m: alignments
        pc0: p(c0) when c0 = ('-', '-')
        pcAll: all p(c), a dict
        b: branch length
        iRate, dRate: insertion deletion rate
    output:
        log of p(m)
    """
    tau = b
    nu = iRate * (tau + 1. / dRate)
    mlen = len(m)
    if mlen == 0:
        logphi = (pc0 - 1) * nu
        logpc = 0
    else:
        logphi = -np.sum(np.log(np.arange(1, mlen+1))) + mlen * np.log(nu) + (pc0 - 1) * nu
        pcvec = [pcAll[c[0]][c[1]] for c in m]
        logpc = np.sum(np.log(pcvec))
    logpm = logphi + logpc
    return logpm


# TODO: EFFICIENCY CAN BE IMPROVED!
def mean_dict_value(pcAllList):
    res = deepcopy(pcAllList[0])
    for ka, va in res.iteritems():
        for kb, vb in pcAllList[0][ka].iteritems():
            res[ka][kb] = np.mean(np.array([pcAllOne[ka][kb] for pcAllOne in pcAllList]))
    return res


# TODO: EFFICIENCY CAN BE IMPROVED!
def logprob_align(b, alignRes, piProb, qMat, iRate, dRate, cList, qRates=[1.]):
    """
    calculate the probability of alignments of segments, or probability of alignments of sequences within one segment
    input:
        alignRes: alignemnt of sequences wthin one segment
        piProb: stationary distribution of qMat
        qMat: rate matrix
        b: branch length
        iRate, dRate: insertion and deletion rate
        cList: list of basic characters
    output:
        logProbAlign: likelihood
    """
    pcAllList = []
    for qRate in qRates:
        qMatExt = q_to_qext(qMat*qRate, dRate)
        pMat = sp.linalg.expm(b*qMatExt)
        piExt = list(piProb) + [0]   # add one element for epsilon
        cListExt = cList + ['-']
        pcAllList.append(prob_c_all(piExt, pMat, b, dRate, cListExt))
    pcAll = mean_dict_value(pcAllList)
    pc0 = pcAll['-']['-']
    logProbAlign = logprob_m(alignRes, pc0, pcAll, b, iRate, dRate)
    return logProbAlign


def felsenstein_peel_recursion(fTildeVec1, fTildeVec2, b1, b2, qMatExt, piProbExt):
    """
    Felsenstain peeling recursion for one step upward
    input:
        fTildeVec1:
        fTildeVec2:
        b1, b2: branch lengths from parent node to two child nodes
        qMatExt: extended version of qMat, with one more column and row for deletion
        piProbExt: extended version of piProb, with one more number 0 for deletino symbol '-'
    output:
        fTildeVec: array, to be used for calculating values of parent node
        fTilde: float, fTildeVec * piProbExt
    """
    pMat1 = sp.linalg.expm(b1 * qMatExt)
    pMat2 = sp.linalg.expm(b2 * qMatExt)
    fTildeVec = np.dot(pMat1, fTildeVec1) * np.dot(pMat2, fTildeVec2)
    fTilde = np.dot(piProbExt, fTildeVec)
    return fTildeVec, fTilde


def ftilde_leaf(piProbExt, ch, nChExt, cListExt):
    """
    fTildeVec, fTilde, beta at leaves: see above for definition
    """
    chIndex = cListExt.index(ch)
    fTildeVec = np.zeros(nChExt)
    fTildeVec[chIndex] = 1
    fTilde = piProbExt[chIndex]
    return fTildeVec, fTilde


def ftilde_recursion_given_tree(msa, seqNames, treeWithBeta, qMatExt, piProbExt, dRate, cListExt):
    """
    probabilisies related to multiple string alignment for a give tree based on PIP for one site
    input:
        msa: multiple alignment, a string of size n, or a list of size n
        seqNames: names of sequences, a list of size n
        treeWithBeta: tree object with .beta method added
        qMatExt, piProbExt, dRate: see above
    output:
        no output, tree is updated with methods .ftilde_vec, .ftilde, and .fv
    """
    nSeq = len(msa)
    nChExt = len(cListExt)
    # isCPhi: indicate if the alignment is all empty.
    isCPhi = msa.count('-') == nSeq
    seqChDict = dict(zip(seqNames, list(msa)))    # sequence name -> character
    # TODO: may further improve recursions for leaves LATER.
    # Set S in PIP paper.
    nonEmptySeqNames = [seqNames[index] for index in xrange(nSeq) if msa[index] != '-']
    # Youngest node in set A in PIP paper.
    if not isCPhi:
        nonZeroFvNodesYoungest = deepcopy(treeWithBeta).mrca(taxon_labels=nonEmptySeqNames)
        node = [node for node in treeWithBeta.postorder_node_iter() if node.index == nonZeroFvNodesYoungest.index][0]
        # Set A in the PIP paper.
        # TODO: replace this with an existing method?
        nonZeroFvNodes = node_path(node, treeWithBeta.seed_node)
    # TODO: improve this part by revising for leaves, internal nodes, and root.
    for node in treeWithBeta.postorder_node_iter():
        if node.is_leaf():
            nodeName = node._get_node_token()
            ch = seqChDict[nodeName]
            fTildeVec, fTilde = ftilde_leaf(piProbExt, ch, nChExt, cListExt)
            beta = node.beta
            if isCPhi:
                fv = 1 + beta * (fTilde - 1)
            else:
                if node in nonZeroFvNodes:
                    fv = beta * fTilde
                else:
                    fv = 0
            # No need to store beta, keep it here for testing purpose.
            node.ftilde_vec = fTildeVec,
            node.ftilde = fTilde
            node.fv = fv
        else:
            child1, child2 = node.child_nodes()
            fTildeVec1 = child1.ftilde_vec[0]
            fTildeVec2 = child2.ftilde_vec[0]
            b1 = child1.edge_length
            b2 = child2.edge_length
            beta = node.beta
            fTildeVec, fTilde = felsenstein_peel_recursion(fTildeVec1, fTildeVec2, b1, b2, qMatExt, piProbExt)
            if isCPhi:
                # Note beta is set to 1 for root.
                fv = 1 + beta * (fTilde - 1)
            else:
                if node in nonZeroFvNodes:
                    # Note beta is set to 1 for root.
                    fv = beta * fTilde
                else:
                    fv = 0
            node.ftilde_vec = fTildeVec,
            node.ftilde = fTilde
            node.fv = fv


def pv_and_beta_given_tree(tree, dRate):
    """
    calculate P(V=v), the prior probability of insertion
    calculate beta(v), survival probability at node given an insertion on the parent edge
    input:
        tree, dRate: see above
    output:
        no output, a new method 'pv' is added to all nodes, with P(V=v) value, and a new method 'beta' is added to all nodes, with beta(v) value.
    """
    tau = tree.length()
    dRateInverse = 1. / dRate
    nuBar = 1. / (tau + dRateInverse)
    iterNodes = tree.preorder_node_iter()
    root = iterNodes.next()    # First one is the root.
    root.pv = dRateInverse * nuBar
    index = 0
    root.index = index
    # For root, set beta = 1, then fv shares the same formula in (4) and (5).
    root.beta = 1
    for node in iterNodes:   # all others
        index += 1
        node.pv = node.edge_length * nuBar
        bMu = node.edge_length * dRate
        if bMu <= 0:
            node.beta = 0
        else:
            node.beta = 1. / bMu * (1 - np.exp(-bMu))
        node.index = index


def pc_from_tree_with_pv_and_fv(treeWithPvAndFv):
    """
    calculate P(C=c), i.e., probability of multiple string alignment of one site using P(V=v) and f_v = P(C=c|V=v)
    input:
        treeWithPvAndFv: tree object with added methods .pv and .fv
    output:
        pc: float, probability
    """
    pc = 0
    for node in treeWithPvAndFv.nodes():
        pc += node.fv * node.pv
    return pc


# TODO: move pv and beta calculation outside to improve efficiency.
def prob_msa_one_site(msa, seqNames, tree, qMatExt, piProbExt, dRate, cListExt):
    """
    calculate P(C=c), i.e., probability of multiple string alignment of one site
    """
    pv_and_beta_given_tree(tree, dRate)
    ftilde_recursion_given_tree(msa, seqNames, tree, qMatExt, piProbExt, dRate, cListExt)
    pc = pc_from_tree_with_pv_and_fv(tree)
    return pc


def logprob_msa_multi_site(msaList, seqNames, tree, qMatExt, piProbExt, dRate, cListExt):
    """
    calculate P(C=c), the alignment log probability of for all loci
    input:
        msaList: alignments, a list of n elements, each of 'ATC' or ('A', 'T', 'C')
        seqNames: name of all sequences, to be consistant with msaList order
        tree:
        qMatExt: extended rate matrix
        piProbExt: extended stationary distribution of qMat
        dRate: deletion rate
        cListExt: extended list of basic characters
    output:
        logProbMsa: array, log probability of multiple string alignments
    """
    nAligns = len(msaList)
    logProbMsa = np.zeros(nAligns)
    index = 0
    for msa in msaList:
        tem = prob_msa_one_site(msa, seqNames, tree, qMatExt, piProbExt, dRate, cListExt)
        logProbMsa[index] = np.log(tem)
        index += 1
    return logProbMsa


def logprob_msa(multiAlign, tree, qMat, piProb, iRate, dRate, cList, qRates=[1.]):
    """
    calculate logP(m) for a multiple alignment under PIP
    1). calculate logP(c) for all unique c in the aligments
    2). construct a dict, c->logP(c), and use this dict for all aligments
    3). calculate logP(m) use the dictionary in 2)
    input:
        multiAlign: dict, seqName->aligned string by loci
        tree: tree, fixed
        qMat, piProb: transition matrix and stationary distribution
        iRate, dRate: insertion and deletion rate
        cList: characher list
    output:
        logProbM: float, log probability of the MSA
    """
    # # Make sure that the tree is rooted.
    # tree.reroot_at_midpoint()
    seqNames = multiAlign.keys()
    nLeaf = len(seqNames)
    msaList = zip(*multiAlign.values())
    mlen = len(msaList)
    msaUnique = list(set(msaList))
    piProbExt = np.append(piProb, 0)
    cListExt = cList + ['-']
    logProbMsaUnique = np.zeros(len(msaUnique))
    for qRate in qRates:
        qMatExt = q_to_qext(qMat * qRate, dRate)
        logProbMsaTem = logprob_msa_multi_site(msaUnique, seqNames, tree, qMatExt, piProbExt, dRate, cListExt)
        logProbMsaUnique += logProbMsaTem
    logProbMsaUnique = logProbMsaUnique / len(qRates)
    logProbMsaDict = dict(zip(msaUnique, logProbMsaUnique))
    logProbMsa = np.array([logProbMsaDict[msa] for msa in msaList])
    # Calculate psi(.,.).
    tau = tree.length()
    nu = iRate * (tau + 1. / dRate)
    cPhi = '-' * nLeaf
    pc0 = prob_msa_one_site(cPhi, seqNames, tree, qMatExt, piProbExt, dRate, cListExt)
    logPsi = -np.sum(np.log(np.arange(1, mlen+1))) + mlen * np.log(nu) + (pc0 - 1) * nu
    logProbM = logPsi + logProbMsa.sum()
    return logProbM


def node_path(nodeDescendent, nodeAncestor):
    """
    node path from nodeDescendent to nodeAncestor
    NOTE: nodeDescendent is included, and nodeAncestor is also included
    NOTE: node.ancestor_iter is equavilent to this function when nodeAncestor is root
    input:
        nodeDescendent: descendent node
        nodeAncestor: ancestror node
    output:
        res: a list of nodes on the path
    """
    res = []
    node = nodeDescendent
    while node is not nodeAncestor:
        res += [node]
        node = node.parent_node
    res += [node]
    return res
