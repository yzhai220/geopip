#!/usr/bin/python

# System modulas.
import sys
import os
import json
from copy import deepcopy
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
from hpip_sim import *
from geopip_sim import *

# Utility functions.
from substitution_util import *
from tree_util import *
from align_util import *

# I/O functions.
from io_json import *
from io_files import *
from io_fasta import *
from io_phylip import *
from io_phyml import *
from io_util import *


cList = ['A', 'C', 'G', 'T']
nList = len(cList)

# Character list: only consider DNA for illustration.
cList = ['A', 'C', 'G', 'T']
nList = len(cList)

# Whether the rate matrix Q is updated or not.
updateQ = True
# Whether a random Q is used in each simulation run.
randomQ = False

# Result folder directory.
resultFolder = dir_main + '/result/hpip_test'
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
rFileLoc = dir_source

###############################################################################
# Start a new simulation using the hPIP model.

# If random Q is used in simulation, then generate a random Q here.
if randomQ:
    qMat = ratem_gtr(4)
    piProb = pi_from_qmat(qMat)
# Otherwise, generate a fixed Q using HKY85 model.
else:
    piProb = np.array([0.25, 0.25, 0.25, 0.25])
    kappa = 2.
    qMat = hky85_from_pi_kappa(piProb, kappa)

# Set deletion rates.
dRate1 = 0.02
dRate2 = 4.
iRateOverdRate = 20

# Generate all indel rates.
ratesList = [(iRateOverdRate*dRate1, dRate1), (iRateOverdRate*dRate2, dRate2)]
nRates = len(ratesList)

# Multinomial probability for two indel rates.
piProbRates = np.array([0.5, 0.5])
# Geometric parameter for generating segments.
p = 0.05

# hPIP papameters
iRateSeg = 2.0
dRateSeg = 0.1
rateSeg = np.array([iRateSeg, dRateSeg])
qMatSeg = np.array([[1., 0], [0, 1.]])
piProbSeg = np.array([0.5, 0.5])

# True tree.
treeStr = '(((seq1:0.05,seq2:0.05):0.05,((seq3:0.05,seq4:0.05):0.05,seq5:0.1):0.05):0.05,((seq6:0.1,seq7:0.05):0.05,seq8:0.25):0.05);'
treeTrue = dendropy.Tree.get_from_string(treeStr, schema="newick")

# Simulate tree from hPIP, for nSeg=20 at root.
nSeg = 20
sim_tree_hpip(treeTrue, iRateSeg, dRateSeg, piProbSeg, qMatSeg, ratesList, piProb, qMat, cList, fixSegNumber=nSeg)
seqs, alignInSeg, alignsSeg, segRateDict, bDict, multiAlignSeqs, lenSegs = leaf_data_from_sim_tree(treeTrue)
multiAlign = multi_align_no_sep(multiAlignSeqs)
print lenSegs
print segRateDict
nSeg = len(segRateDict.keys())

# Unroot the true tree and then output the true tree.
treeTrue.deroot()
write_tree(treeTrue, dataLoc + '/tree_true.txt')

# Output the fasta file.
dict_write_fasta(multiAlign, dataLoc)
# Output the phylip file.
dict_write_phylip(multiAlign, dataLoc)
# Output the phylip file with only substitutions.
dict_write_phylip_subs_only(multiAlign, dataLoc)

# Output true values used in simulation.
write_dist_from_bdict(bDict, dataLoc)
output_lensegs(lenSegs, dataLoc)
output_ratelist(ratesList, dataLoc)
output_qmat(qMat, dataLoc)

###############################################################################
# Start the inference steps.

# phyml.
phymlLoc = dir_main + '/software/phyml'

time1 = time()
# Q from phyml is used as start value in other methods.
phyml_run(dataLoc, phymlLoc, fileName='msa.phylip.txt')
time1 = time() - time1
qMatPhyml = get_qmat_from_phyml_output(dataLoc, fileName='stats_phyml_1rate.txt')
qMatStart = qMatPhyml

# CTMC+NJ.
time2 = time()
qMatCTMC1, bDictCTMC1, treeCTMC1, nllkCTMC1 = opt_ctmc_full(qMatStart, multiAlign, javaDirectory, modelDirectory, eStepFile, parametersPath, inputLoc, outputLoc, dataLoc, execsLoc, rFileLoc, cList, qRates=[1.], suffix='_ctmc_1rate', updateQ=updateQ, tol=1.e-2, bTol=1.e-3, iterMax=100)
time2 = time() - time2

# Output results.
write_dist_from_bdict(bDictCTMC1, dataLoc, suffix='_ctmc_1rate')
output_qmat(qMatCTMC1, dataLoc, suffix='_ctmc_1rate')

# PIP+NJ.
# Random start rates.
rateStart = pip_start_data(multiAlign)
time3 = time()
rateEstPIP, qMatEstPIP, bDictEstPIP, treeEstPIP = opt_pip_full(rateStart, qMatStart, multiAlign, javaDirectory, modelDirectory, eStepFile, parametersPath, inputLoc, outputLoc, dataLoc, execsLoc, rFileLoc, cList, qRates=[1.], suffix='_pip', updateQ=updateQ, updateRate=True, updateRateFixdRateTimesb=True, tol=1.e-2, bTol=1.e-3, iterMax=100)
time3 = time() - time3

# Output results.
write_dist_from_bdict(bDictEstPIP, dataLoc, suffix='_pip')
output_qmat(qMatEstPIP, dataLoc, suffix='_pip')
output_ratelist(rateEstPIP, dataLoc, suffix='_pip')

# GeoPIP3+NJ.
m = 3
# Arbitrary starts.
pStart = 0.1
piProbRatesStart = np.ones(m) / m
iRateOverdRate = 20
lenSegsStart = [len(multiAlign.values()[0])]

# Random starts.
ratesListStart = []
for i in xrange(m):
    dRate = rd.random() + i
    iRate = dRate * iRateOverdRate
    ratesListStart.append((iRate, dRate))

dRateMean = sum(zip(*ratesListStart)[1]) / m
iRateMean = dRateMean * iRateOverdRate
segRateDictStart = {0: (iRateMean, dRateMean)}

suffix1 = '_geopip_' + str(m)

time4 = time()
bDictEst1, qMatEst1, segRateDictEst1, alignsInSegEst1, pEst1, piProbRatesEst1, ratesListEst1, lenSegsEst1, treeEst1, nllk1 = opt_geopip_full(m, pStart, qMatStart, segRateDictStart, piProbRatesStart, ratesListStart, multiAlign, lenSegsStart, javaDirectory, modelDirectory, eStepFile, parametersPath, inputLoc, outputLoc, dataLoc, execsLoc, rFileLoc, cList, suffix=suffix1, updateQ=updateQ, updateSeg=True, updateRate=True, updateRateFixdRateTimesTau=True, tol=1.e-2, bTol=1.e-3, iterMax=100)
time4 = time() - time4

# Output results.
write_dist_from_bdict(bDictEst1, dataLoc, suffix=suffix1)
output_qmat(qMatEst1, dataLoc, suffix=suffix1)
output_lensegs(lenSegsEst1, dataLoc, suffix=suffix1)
output_ratelist(ratesListEst1, dataLoc, suffix=suffix1)

rateCatListEst1 = ratecatlist_from_segratedict_and_ratelist(segRateDictEst1, ratesListEst1)
output_ratecat(rateCatListEst1, dataLoc, suffix=suffix1)
rateCatAtEachLoc = output_ratecatateachloc(segRateDictEst1, ratesListEst1, lenSegsEst1, dataLoc, suffix=suffix1)

output_rate_with_msa(multiAlign, rateCatAtEachLoc, dataLoc + '/msa_rate' + suffix1 + '.txt')


# GeoPIP5+NJ.
m = 5
# Arbitrary starts.
pStart = 0.1
piProbRatesStart = np.ones(m) / m
iRateOverdRate = 20
lenSegsStart = [len(multiAlign.values()[0])]
# Random starts.
ratesListStart = []
for i in xrange(m):
    dRate = rd.random() + i
    iRate = dRate * iRateOverdRate
    ratesListStart.append((iRate, dRate))

dRateMean = sum(zip(*ratesListStart)[1]) / m
iRateMean = dRateMean * iRateOverdRate
segRateDictStart = {0: (iRateMean, dRateMean)}

suffix2 = '_geopip_' + str(m)

time5 = time()
bDictEst2, qMatEst2, segRateDictEst2, alignsInSegEst2, pEst2, piProbRatesEst2, ratesListEst2, lenSegsEst2, treeEst2, nllk2 = opt_geopip_full(m, pStart, qMatStart, segRateDictStart, piProbRatesStart, ratesListStart, multiAlign, lenSegsStart, javaDirectory, modelDirectory, eStepFile, parametersPath, inputLoc, outputLoc, dataLoc, execsLoc, rFileLoc, cList, suffix=suffix2, updateQ=updateQ, updateSeg=True, updateRate=True, updateRateFixdRateTimesTau=True, tol=1.e-2, bTol=1.e-3, iterMax=100)
time5 = time() - time5

# Output results.
write_dist_from_bdict(bDictEst2, dataLoc, suffix=suffix2)
output_qmat(qMatEst2, dataLoc, suffix=suffix2)
output_lensegs(lenSegsEst2, dataLoc, suffix=suffix2)
output_ratelist(ratesListEst2, dataLoc, suffix=suffix2)

rateCatListEst2 = ratecatlist_from_segratedict_and_ratelist(segRateDictEst2, ratesListEst2)
output_ratecat(rateCatListEst2, dataLoc, suffix=suffix2)
rateCatAtEachLoc = output_ratecatateachloc(segRateDictEst2, ratesListEst2, lenSegsEst2, dataLoc, suffix=suffix2)

output_rate_with_msa(multiAlign, rateCatAtEachLoc, dataLoc + '/msa_rate' + suffix2 + '.txt')


###############################################################################
# Post-processing.

# All trees considered.
treeNames = ['true', 'ctmc_1rate', 'phyml_1rate', 'pip', suffix1[1:], suffix2[1:]]
tns = dendropy.TaxonNamespace()

# Get all trees.
treeAllNonScaled = []
treeAllScaled = []
for treeName in treeNames:
    treeFile = dataLoc+'/tree_' + treeName + '.txt'
    treeTemp = dendropy.Tree.get_from_path(treeFile, schema='newick', taxon_namespace=tns)
    treeAllNonScaled.append(treeTemp)
    treeTemp = dendropy.Tree.get_from_path(treeFile, schema='newick', taxon_namespace=tns)
    treeTemp.scale_edges(1. / treeTemp.length())
    treeAllScaled.append(treeTemp)

# Output wRF among unscaled trees.
treeDictNonScaled = dict(zip(treeNames, treeAllNonScaled))
treeDistNonScaled = dist_among_trees(treeDictNonScaled)
outDistNonScaledFile = dataLoc + '/dist_trees_non_scaled.txt'
with open(outDistNonScaledFile, 'w') as f:
    json.dump(treeDistNonScaled, f, sort_keys=True, indent=4)

# Output wRF among scaled trees.
treeDictScaled = dict(zip(treeNames, treeAllScaled))
treeDistScaled = dist_among_trees(treeDictScaled)
outDistScaledFile = dataLoc + '/dist_trees_scaled.txt'
with open(outDistScaledFile, 'w') as f:
    json.dump(treeDistScaled, f, sort_keys=True, indent=4)

# Output RF among trees.
treeDistSym = dist_among_trees_sym(treeDictNonScaled)
outDistSymFile = dataLoc + '/dist_trees_sym.txt'
with open(outDistSymFile, 'w') as f:
    json.dump(treeDistSym, f, sort_keys=True, indent=4)


####################################################
# Output result summary into file print.txt.
outputFile = dataLoc + '/print.txt'
sys.stdout = open(outputFile, 'w')

print '### summary statistics ###'
print 'number of sequences:', len(multiAlign)
print 'number of sites:', len(multiAlign.values()[0])
print 'nSeg:', nSeg
print 'indel rates:', ratesList
print 'hPIP rates', iRateSeg, dRateSeg
print 'stationary distribution:', piProb
print 'Q:', qMat
print 'time order: phyml, ctmc, pip, geopip3, geopip5'
print 'time:', time1, time2, time3, time4, time5

print '### Comparison of pariwise branch lengths ###'
print 'pair, true tree, CTMC tree, PIP tree, GeoPIP tree'
for key in bDict.keys():
    print key, bDict[key], bDictCTMC1[key], bDictEstPIP[key], bDictEst1[key], bDictEst2[key]

print '### Weighted Robinson Foulds (wRF) of unscaled trees and true tree ###'
print json.dumps(treeDistNonScaled['true'], sort_keys=True, indent=4)

print '### Weighted Robinson Foulds (wRF) of scaled trees and true tree ###'
print json.dumps(treeDistScaled['true'], sort_keys=True, indent=4)

print '### Robinson Foulds (RF) of scaled trees trees and true tree ###'
print json.dumps(treeDistSym['true'], sort_keys=True, indent=4)
