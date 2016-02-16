# Functions for operations on the substitution rate matrix.
# This is based on the weight and feature generation of the rate matrix.

import numpy as np
from itertools import combinations
from copy import deepcopy


def f81_from_pi(piProb):
    """
    generate rate matrix qMat from stationary dist. piProb, using F81 model
    rate matrix is scaled
    """
    piProb = np.array(piProb)
    n = len(piProb)
    res = np.tile(piProb, (n, 1))
    np.fill_diagonal(res, 0)
    rowsum = res.sum(axis=1)
    np.fill_diagonal(res, -rowsum)
    beta = 1.0 / (1. - sum(piProb**2))
    res = beta * res
    return res


def hky85_from_pi_kappa(piProb, kappa):
    """
    generate rate matrix qMat from stationary dist. piProb, and parameter kappa, using HKY85 model
    """
    piProb = np.array(piProb)
    n = len(piProb)
    res = np.tile(piProb, (n, 1))
    res[0, 2] = res[0, 2] * kappa   # assuming 4*4 Q
    res[1, 3] = res[1, 3] * kappa   # assuming 4*4 Q
    res[2, 0] = res[2, 0] * kappa   # assuming 4*4 Q
    res[3, 1] = res[3, 1] * kappa   # assuming 4*4 Q
    np.fill_diagonal(res, 0)
    rowsum = res.sum(axis=1)
    np.fill_diagonal(res, -rowsum)
    beta = 1.0 / (2 * (piProb[0] + piProb[2]) * (piProb[1] + piProb[3]) + 2 * kappa * (piProb[0] * piProb[2] + piProb[1] * piProb[3]))
    res = beta * res
    return res


def unif_gtr(nCharType):
    """
    univariate features for GTR model
    """
    wLen = (nCharType + 1) * nCharType / 2
    uniF = np.zeros((nCharType, wLen))
    np.fill_diagonal(uniF, 1)
    return uniF


def bivariatef_gtr(nCharType):
    """
    generate bivariate features for GTR model
    """
    wLen = (nCharType + 1) * nCharType / 2
    bivariateF = np.zeros((wLen, nCharType, nCharType))
    biIters = combinations(range(nCharType), 2)
    for i in range(nCharType, wLen):
        j, k = biIters.next()
        bivariateF[i, j, k] = 1
        bivariateF[i, k, j] = 1
    return bivariateF


def feature_gtr(nCharType):
    """
    generate feature matrices for GTR model
    """
    feature = bivariatef_gtr(nCharType)
    for i in range(nCharType):
        feature[i, :, i] = 1.
    return feature


def ratem_from_w_grad(w, grad):
    """
    calculate rate matrix qMat from weight w and gradient grad
    """
    qMat = deepcopy(grad)
    wLen = w.size
    for i in range(wLen):
        qMat[i, :, :] = grad[i, :, :] * w[i]
    qMat = qMat.sum(axis=0)
    qMat = np.exp(qMat)
    np.fill_diagonal(qMat, 0)
    np.fill_diagonal(qMat, -qMat.sum(axis=1))
    return qMat


def pi_from_w_unif(w, uniF):
    """
    stationary distribution from weight and univariate features
    """
    piProbUnscaled = np.exp(uniF.dot(w))
    piSum = piProbUnscaled.sum()
    piScaled = piProbUnscaled / piSum
    return piScaled


def ratem_norm_from_w_unif_grad(w, uniF, grad):
    """
    normalized rate matrix, from w, uniF, gradient
    """
    piProb = pi_from_w_unif(w, uniF)
    qMat = ratem_from_w_grad(w, grad)
    beta = -1. / (piProb * qMat.diagonal()).sum()
    qMatNorm = beta * qMat
    return qMatNorm


def ratem_gtr(nCharType=4):
    """
    generate a rate matrix Q for a given dimension
    """
    nWeight = nCharType * (nCharType + 1) / 2
    w = np.random.normal(size=nWeight)
    grad = feature_gtr(nCharType)
    uniF = unif_gtr(nCharType)
    qMat = ratem_norm_from_w_unif_grad(w, uniF, grad)
    return qMat
