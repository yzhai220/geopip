#!/usr/bin/python

# System modulas.
import sys
import os
from time import time

# Get directories.
dir_script = os.getcwd()
dir_main = dir_script[:dir_script.rfind('/')]
dir_source = dir_main + '/src'
# Add source path into import path.
sys.path.append(dir_source)

# Estimation functions.
from pip_est import *
from ctmc_est import *
from geopip_est import *
from geopip_sim import *

# Utility functions.
from substitution_util import *
from tree_util import *
from align_util import *
from seg_util import *

# I/O functions.
from io_json import *
from io_files import *
from io_fasta import *
from io_phylip import *
from io_phyml import *
from io_util import *

# Character list: only consider DNA for illustration.
cList = ['A', 'C', 'G', 'T']
nList = len(cList)

# Result folder directory.
resultFolder = dir_main + '/result/segmentation_test'
if not os.path.exists(resultFolder):
    os.makedirs(resultFolder)
os.chdir(resultFolder)
# Substitution model directory.
modelDirectory = dir_main + '/model/gtr-model'
# Result folders for estimating rate matrix using Java code.
eStepFile, parametersPath, inputLoc, outputLoc, dataLoc = make_global_location_files(resultFolder, 'runs')
javaDirectory = dir_main + '/software/java'
# Result folder for dumping Java outputs.
execsLoc = resultFolder + '/state/execs'
if not os.path.exists(execsLoc):
    os.makedirs(execsLoc)
# Exact result folder for every run.
subFolderPath = resultFolder + '/runs'
if not os.path.exists(subFolderPath):
    os.makedirs(subFolderPath)

###############################################################################
# Start a new simulation for segmentation analysis.
# Generate a random Q rate matrix using GTR model.
qMat = ratem_gtr(nList)
piProb = pi_from_qmat(qMat)

# Set deletion rates.
dRate1 = 0.02
dRate2 = 2.
# Set insertion rates by deletion rate * iRateOverdRate.
iRateOverdRate = 20
ratesList = [(iRateOverdRate*dRate1, dRate1), (iRateOverdRate*dRate2, dRate2)]
nRates = len(ratesList)
# Set multinomial probability for rate categories.
piProbRates = np.ones(len(ratesList)) / float(sum(np.ones(len(ratesList))))

# Set geometric parameter.
p = 0.05
# Set number of leaves, 2**log2Loeaves.
log2Leaves = 4

# Generate a perfect binary tree, with 2**log2Loeaves and edge length bLen.
treeStr = sim_tree_newick_fixed_edge_length(log2Leaves, bLen=0.1)
treeTrue = dendropy.Tree.get_from_string(treeStr, schema="newick")

# Number of segments generated at the root.
nSeg = 20
# Simulate data using the GeoPIP model.
sim_tree(treeTrue, p, ratesList, piProbRates, piProb, qMat, cList, fixSegNumber=nSeg)
# Keep only data at leaves for inference.
seqs, alignInSeg, alignsSeg, segRateDict, bDict, multiAlignSeqs, lenSegs = leaf_data_from_sim_tree(treeTrue)
# Get multiple sequence alignment.
multiAlign = multi_align_no_sep(multiAlignSeqs)

# Record running times.
timeVec = []

# Output results at print.txt.
outputFile = dataLoc + '/print.txt'
sys.stdout = open(outputFile, 'w')

print 'error rate, sequence used:'

# Record true indel rates for all loci, for comparison.
indelRatesAll = indel_rate_for_all_loci(segRateDict, ratesList, lenSegs)

# Number of sequences used for segmetnation, n = 2**log2n.
log2n = 1
# Get the subtree with only sequences picked for segmentation.
suffixPre = '_sub_' + str(2**log2n)
seqNames = ['seq' + str(i) for i in range(1, 2**log2n+1)]
treeSub = subtree_retain_taxa_with_labels(treeTrue, seqNames)

# Get the multiple sequence alignment for only selected sequences.
multiAlignSub, lenSegsSub, msaColIndexSmall = msa_retain_taxa_with_labels_and_track_col_num(multiAlign, seqNames, lenSegs)
msaListSub = zip(*multiAlignSub.values())

# Running the segmentation algorithm.
timeStart = time()
lenSegsSubEst, rateSegsSubEst = mle_seg_len_rate(p, msaListSub, ratesList, seqNames, treeSub, qMat, piProb, piProbRates, cList)
timeEnd = time()
timeVec.append(timeEnd - timeStart)

# Record estimated rates.
msaColIndex = msaColIndexSmall
# Print results.
print prop_of_col_with_wrong_rates(indelRatesAll, msaColIndexSmall, msaColIndex, lenSegsSubEst, rateSegsSubEst), '_to_'.join([seqNames[0], seqNames[-1]])

# Use more sequences, by a facter of 2, until all sequences are used.
for i in xrange(log2Leaves-1):
    log2n = i + 2
    suffixPre = '_sub_' + str(2**log2n)
    seqNames = ['seq' + str(i) for i in range(1, 2**log2n+1)]
    treeSub = subtree_retain_taxa_with_labels(treeTrue, seqNames)
    multiAlignSub, lenSegsSub, msaColIndex = msa_retain_taxa_with_labels_and_track_col_num(multiAlign, seqNames, lenSegs)
    msaListSub = zip(*multiAlignSub.values())
    timeStart = time()
    lenSegsSubEst, rateSegsSubEst = mle_seg_len_rate(p, msaListSub, ratesList, seqNames, treeSub, qMat, piProb, piProbRates, cList)
    timeEnd = time()
    timeVec.append(timeEnd - timeStart)
    print prop_of_col_with_wrong_rates(indelRatesAll, msaColIndexSmall, msaColIndex, lenSegsSubEst, rateSegsSubEst), '_to_'.join([seqNames[0], seqNames[-1]])

print 'indel rates:'
print ratesList
print 'Segmentation time:'
print timeVec
