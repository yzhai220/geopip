# Input and output functions for using PhyML.

import os
import numpy as np
import scipy as sp
from scipy.stats import chi2
from subprocess import Popen, PIPE
from shlex import split

from pip_util import pi_from_qmat


# TODO: THIS NEED TO BE CHANGED.
# phymlPre1 = 'phyml -d nt -m GTR -c 1'
# phymlPre4 = 'phyml -d nt -m GTR'


def phyml_run(dataLoc, phymlLoc, fileName='msa.phylip.txt', suffix=''):
    """
    use Phyml software to estimate tree
    input:
        dataLoc: data folder
        fileName: alignment file, by default 'msa.phylip.txt', in phylip format
    output:
        'tree_phyml_1rate.txt', 'stats_phyml_1rate.txt': no rate variation
        'tree_phyml_4rate.txt', 'stats_phyml_4rate.txt': rate variation with 4
            rates
    """
    phymlPre1 = phymlLoc + ' -d nt -m GTR -c 1'
    phymlPre4 = phymlLoc + ' -d nt -m GTR'
    alignFile = dataLoc + '/' + fileName
    phymlStr1 = phymlPre1 + ' -i ' + alignFile
    phymlStr4 = phymlPre4 + ' -i ' + alignFile
    outputStats = dataLoc + '/' + fileName + '_phyml_stats.txt'
    outputTree = dataLoc + '/' + fileName + '_phyml_tree.txt'
    os.system(phymlStr1)
    outputStatsRename = dataLoc + '/' + 'stats_phyml_1rate' + suffix + '.txt'
    statsRenameStr = 'mv ' + outputStats + ' ' + outputStatsRename
    outputTreeRename = dataLoc + '/' + 'tree_phyml_1rate' + suffix + '.txt'
    treeRenameStr = 'mv ' + outputTree + ' ' + outputTreeRename
    os.system(statsRenameStr)
    os.system(treeRenameStr)
    os.system(phymlStr4)
    outputStatsRename = dataLoc + '/' + 'stats_phyml_4rate' + suffix + '.txt'
    statsRenameStr = 'mv ' + outputStats + ' ' + outputStatsRename
    outputTreeRename = dataLoc + '/' + 'tree_phyml_4rate' + suffix + '.txt'
    treeRenameStr = 'mv ' + outputTree + ' ' + outputTreeRename
    os.system(statsRenameStr)
    os.system(treeRenameStr)


def get_qmat_from_phyml_output(dataLoc, fileName='stats_phyml_4rate.txt'):
    outputFile = dataLoc + '/' + fileName
    f = open(outputFile)
    lines = f.readlines()
    f.close()
    qMatIndexStart = lines.index('  [A---------C---------G---------T------]\n') + 1
    qMatIndexEnd = qMatIndexStart + 4
    qMatLines = lines[qMatIndexStart:qMatIndexEnd]
    qMat = []
    for line in qMatLines:
        tem = line.split()
        temNum = [float(item) for item in tem]
        qMat.append(temNum)
    qMat = np.array(qMat)
    return qMat


def get_alpha_from_phyml_output(dataLoc, fileName='stats_phyml_4rate.txt'):
    outputFile = dataLoc + '/' + fileName
    f = open(outputFile)
    lines = f.readlines()
    f.close()
    res = [line for line in lines if line.startswith('  - Gamma shape parameter:')]
    alpha = float(res[0].split()[-1])
    return alpha


# assuming 4 * 4
def qmat_paml_format(qMat):
    """
    get PaML format from a rate matrix, i.e., frequencies and relative rate parameters
    input:
        qMat: a rate matrix
    output:
        piProb: stationary distribution of the rate matrix
        rate: relative rate parameters, in the order of
              A <-> C, A <-> G, A <-> T, C <-> G, C <-> T, G <-> T
    """
    piProb = pi_from_qmat(qMat)
    tMat = qMat / piProb
    # dimRow, dimCol = qMat.shape
    # tem = np.array([tMat[i, (i+1):4] for i in xrange(4-1)])
    rate = list([list(tMat[0, 1:4]), list(tMat[1, 2:4]), list(tMat[2, 3:4])])
    rate = np.hstack(rate)
    rate = rate / rate[-1]
    return piProb, rate


def write_phyml_input(qMat, dataLoc, rate1=False, alignFileName='msa.phylip.txt', inputFileName='phyml_input.txt'):
    """
    wriet a txt file containing input information for running Phyml using Phylip style interface and pipes
    """
    alignFile = dataLoc + '/' + alignFileName
    inputFile = dataLoc + '/' + inputFileName
    piProb, rate = qmat_paml_format(qMat)
    f = open(inputFile, 'w')
    f.write(alignFile+'\n')
    f.write('+\n')
    [f.write('M\n') for i in xrange(4)]    # change to custom model
    f.write('E\n')    # stationary distribution given
    [f.write(piProb[i].__str__()+'\n') for i in xrange(len(piProb))]
    f.write('\n')    # this confirms regularization of pi to 1
    f.write('K\n')    # set custom model relative changing ratio
    f.write('012345\n')    # set the custom model
    [f.write(rate[i].__str__()+'\n') for i in xrange(len(rate))]   # set relative ratio
    f.write('O\n')    # fix relative changing ratio
    if rate1:
        f.write('R\n')    # only one rate categery
    f.write('Y\n')    # fix relative changing ratio
    f.close()


def phyml_run_fix_q(dataLoc, rate1=False, alignFileName='msa.phylip.txt', inputFileName='phyml_input.txt', outputFileName=''):
    """
    use Phyml software to estimate tree, when Q is fixed
    """
    inputFile = dataLoc + '/' + inputFileName
    phymlStr = 'cat ' + inputFile
    p1 = Popen(split(phymlStr), stdout=PIPE)
    p2 = Popen(split("phyml"), stdin=p1.stdout)
    p2.wait()
    outputStats = dataLoc + '/' + alignFileName + '_phyml_stats.txt'
    outputTree = dataLoc + '/' + alignFileName + '_phyml_tree.txt'
    if rate1:
        outputStatsRename = dataLoc + '/' + 'stats_phyml_1rate' + outputFileName + '.txt'
        outputTreeRename = dataLoc + '/' + 'tree_phyml_1rate' + outputFileName + '.txt'
    else:
        outputStatsRename = dataLoc + '/' + 'stats_phyml_4rate' + outputFileName + '.txt'
        outputTreeRename = dataLoc + '/' + 'tree_phyml_4rate' + outputFileName + '.txt'
    statsRenameStr = 'mv ' + outputStats + ' ' + outputStatsRename
    treeRenameStr = 'mv ' + outputTree + ' ' + outputTreeRename
    os.system(statsRenameStr)
    os.system(treeRenameStr)


def write_phyml_input_fix_alpha(qMat, dataLoc, alpha, alignFileName='msa.phylip.txt', inputFileName='phyml_input.txt'):
    """
    wriet a txt file containing input information for running Phyml using Phylip style interface and pipes
    """
    alignFile = dataLoc + '/' + alignFileName
    inputFile = dataLoc + '/' + inputFileName
    piProb, rate = qmat_paml_format(qMat)
    f = open(inputFile, 'w')
    f.write(alignFile+'\n')
    f.write('+\n')
    [f.write('M\n') for i in xrange(4)]    # change to custom model
    f.write('E\n')    # stationary distribution given
    [f.write(piProb[i].__str__()+'\n') for i in xrange(len(piProb))]
    f.write('\n')    # this confirms regularization of pi to 1
    f.write('K\n')    # set custom model relative changing ratio
    f.write('012345\n')    # set the custom model
    [f.write(rate[i].__str__()+'\n') for i in xrange(len(rate))]   # set relative ratio
    f.write('O\n')    # fix relative changing ratio
    f.write('A\n')
    f.write('N\n')
    f.write(str(alpha)+'\n')
    f.write('Y\n')    # fix relative changing ratio
    f.close()


def cutting_point(k, alpha=0.5):
    """
    calculate cutting points for k categories of Gamma(alpha, beta=alpha)
    output: an array of k-1 cutting values, excluding 0 and inf
    """
    a = np.array(range(1, k), dtype=float) / k
    res = chi2.ppf(a, 2*alpha) / (2*alpha)
    return res


def mean_values(k, alpha=0.5):
    """
    calculate mean values for k categories of Gamma(alpha, beta=alpha)
    output: an array of k mean values
    """
    cuttingPoints = cutting_point(k, alpha)
    cons = 1. / k
    temVec = cuttingPoints * alpha
    intVec = sp.special.gammainc(alpha+1, temVec)
    intVec = np.append(np.append(0, intVec), 1)
    r = np.diff(intVec) / cons
    r = list(r)
    return r


def median_values(k, alpha=0.5):
    """
    calculate median values for k categories of Gamma(alpha, beta=alpha)
    output: an array of k mean values
    """
    cuttingPoints = cutting_point(2*k, alpha)
    index = np.array(range(k)) * 2
    r = cuttingPoints[index]
    r = r / r.mean()
    return r
