# Basic operation functions used in PIP.

import numpy as np


def pi_from_qmat(qMat):
    """
    Calculate stationary distribution pi from qMat
    Qext * x = b where b = [0, 0, ..., 1]
    Qext = Q add one row with all 1s
    input:
        qMat: a transition matrix
    output:
        x: stationary distritbution of qMat
    """
    qMatDim = len(qMat)
    if qMatDim == 1:
        x = np.array([1.])
    else:
        qMatExt = np.hstack([qMat, np.ones((qMatDim, 1.0))])
        b = np.zeros((qMatDim+1, ))
        b[-1] = 1
        res = np.linalg.lstsq(qMatExt.T, b)
        x = res[0]
        # x = list(x)
    return x


def q_to_qext(qMat, dRate):
    """
    Transform Q to Q* with one more column and row representing deletion
    input:
        qMat: rate matrix
        dRate: deletion rate
    output:
        qMatExt: extended version of qMat
    """
    aShape = np.shape(qMat)
    n = aShape[0]
    qMatExt = np.zeros((n+1, n+1), dtype=float)
    qMatExt[:n, :n] = qMat
    qMatExt[:n, n] = dRate
    np.fill_diagonal(qMatExt, 0)
    rowSum = qMatExt.sum(axis=1)
    np.fill_diagonal(qMatExt, -rowSum)
    return qMatExt
