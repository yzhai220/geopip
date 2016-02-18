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

# Character list: only consider DNA for illustration.
cList = ['A', 'C', 'G', 'T']
nList = len(cList)

# Whether the rate matrix Q is updated or not.
updateQ = True

# Result folder directory.
resultFolder = dir_main + '/result/geopip_test'
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
# Start a new simulation using the GeoPIP model.

# Generate a random Q using GTR model
qMat = ratem_gtr(nList)
piProb = pi_from_qmat(qMat)

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

# Generate a perfect binary tree with 2**log2nleaf leaves and edge length bLen.
log2nleaf = 3   # 2**3 = 8 leaves.
bLen = 0.1
treeStr = sim_tree_newick_fixed_edge_length(log2nleaf, bLen)
treeTrue = dendropy.Tree.get_from_string(treeStr, schema="newick")

# Simulation using the GeoPIP model.
# Number of segments, to be fixed as nSeg.
nSeg = 5
# A geometric number of segments will be generated if nSeg is not provided,
sim_tree(treeTrue, p, ratesList, piProbRates, piProb, qMat, cList, fixSegNumber=nSeg)
# Keep only sequences at leaves for inference.
seqs, alignInSeg, alignsSeg, segRateDict, bDict, multiAlignSeqs, lenSegs = leaf_data_from_sim_tree(treeTrue)
# Get the true MSA after removing empty columns.
multiAlign = multi_align_no_sep(multiAlignSeqs)

# Number of segments generated.
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

# 1. phyml inference.
# Record running time.
# phyml
phymlLoc = dir_main + '/software/phyml'

time1 = time()
phyml_run(dataLoc, phymlLoc, fileName='msa.phylip.txt')
time1 = time() - time1

# Use Q from phyml as starting value qMatStart
qMatPhyml = get_qmat_from_phyml_output(dataLoc, fileName='stats_phyml_1rate.txt')
alpha = get_alpha_from_phyml_output(dataLoc)
qMatStart = qMatPhyml

# 2. CTMC+NJ.
time2 = time()
qMatCTMC1, bDictCTMC1, treeCTMC1, nllkCTMC1 = opt_ctmc_full(qMatStart, multiAlign, javaDirectory, modelDirectory, eStepFile, parametersPath, inputLoc, outputLoc, dataLoc, execsLoc, rFileLoc, cList, qRates=[1.], suffix='_ctmc_1rate', updateQ=updateQ, tol=1.e-2, bTol=1.e-3, iterMax=100)
time2 = time() - time2

# Output results.
write_dist_from_bdict(bDictCTMC1, dataLoc, suffix='_ctmc_1rate')
output_qmat(qMatCTMC1, dataLoc, suffix='_ctmc_1rate')

# 3. PIP+NJ
# Random start rates.
rateStart = pip_start_data(multiAlign)
time3 = time()
rateEstPIP, qMatEstPIP, bDictEstPIP, treeEstPIP = opt_pip_full(rateStart, qMatStart, multiAlign, javaDirectory, modelDirectory, eStepFile, parametersPath, inputLoc, outputLoc, dataLoc, execsLoc, rFileLoc, cList, qRates=[1.], suffix='_pip', updateQ=updateQ, updateRate=True, updateRateFixdRateTimesb=True, tol=1.e-2, bTol=1.e-3, iterMax=100)
time3 = time() - time3

# Output results.
write_dist_from_bdict(bDictEstPIP, dataLoc, suffix='_pip')
output_qmat(qMatEstPIP, dataLoc, suffix='_pip')
output_ratelist(rateEstPIP, dataLoc, suffix='_pip')

# 4. GeoPIP+NJ with true segmentation start.
# The number of rate categories m is fixed as 2.
m = 2
# Start values for the geometric parameters p, and multinomial probability pi.
pStart = 0.1
piProbRatesStart = np.array([0.5, 0.5])

# Suffix for output results.
geopipSuffix = '_geopip_true_start'

# Start values: true segmentation start.
segRateDictStart = segRateDict
ratesListStart = ratesList
lenSegsStart = lenSegs

# Run GeoPIP+NJ.
time4 = time()
bDictEst1, qMatEst1, segRateDictEst1, alignsInSegEst1, pEst1, piProbRatesEst1, ratesListEst1, lenSegsEst1, treeEst1, nllk1 = opt_geopip_full(m, pStart, qMatStart, segRateDictStart, piProbRatesStart, ratesListStart, multiAlign, lenSegsStart, javaDirectory, modelDirectory, eStepFile, parametersPath, inputLoc, outputLoc, dataLoc, execsLoc, rFileLoc, cList, suffix=geopipSuffix, updateQ=updateQ, updateSeg=True, updateRate=True, updateRateFixdRateTimesTau=True, tol=1.e-2, bTol=1.e-3, iterMax=100)
time4 = time() - time4

# Output results.
write_dist_from_bdict(bDictEst1, dataLoc, suffix=geopipSuffix)
output_qmat(qMatEst1, dataLoc, suffix=geopipSuffix)
output_lensegs(lenSegsEst1, dataLoc, suffix=geopipSuffix)
output_ratelist(ratesListEst1, dataLoc, suffix=geopipSuffix)

rateCatListEst1 = ratecatlist_from_segratedict_and_ratelist(segRateDictEst1, ratesListEst1)
output_ratecat(rateCatListEst1, dataLoc, suffix=geopipSuffix)
rateCatAtEachLoc = output_ratecatateachloc(segRateDictEst1, ratesListEst1, lenSegsEst1, dataLoc, suffix=geopipSuffix)
# Output estimated rate category indicators together with the MSA.
output_rate_with_msa(multiAlign, rateCatAtEachLoc, dataLoc + '/msa_rate_geopip_true_start.txt')

# 5. GeoPIP+NJ with random start.
# Random start values.
iRateOverdRate = 20
lenSegsStart = [len(multiAlign.values()[0])]
dRate1 = rd.random()
iRate1 = dRate1 * iRateOverdRate
dRate2 = rd.random() + 1.
iRate2 = dRate2 * iRateOverdRate
ratesListStart = [(iRate1, dRate1), (iRate2, dRate2)]

dRateMean = (dRate1 + dRate2) / 2
iRateMean = dRateMean * iRateOverdRate
# Start segmentation with only one segment for all MSA columns.
segRateDictStart = {0: (iRateMean, dRateMean)}

# Run GeoPIP+NJ>
time5 = time()
bDictEst2, qMatEst2, segRateDictEst2, alignsInSegEst2, pEst2, piProbRatesEst2, ratesListEst2, lenSegsEst2, treeEst2, nllk2 = opt_geopip_full(m, pStart, qMatStart, segRateDictStart, piProbRatesStart, ratesListStart, multiAlign, lenSegsStart, javaDirectory, modelDirectory, eStepFile, parametersPath, inputLoc, outputLoc, dataLoc, execsLoc, rFileLoc, cList, suffix='_geopip', updateQ=updateQ, updateSeg=True, updateRate=True, updateRateFixdRateTimesTau=True, tol=1.e-2, bTol=1.e-3, iterMax=100)
time5 = time() - time5

# Output results.
write_dist_from_bdict(bDictEst2, dataLoc, suffix='_geopip')
output_qmat(qMatEst2, dataLoc, suffix='_geopip')
output_lensegs(lenSegsEst2, dataLoc, suffix='_geopip')
output_ratelist(ratesListEst2, dataLoc, suffix='_geopip')

rateCatListEst2 = ratecatlist_from_segratedict_and_ratelist(segRateDictEst2, ratesListEst2)
output_ratecat(rateCatListEst2, dataLoc, suffix='_geopip')
rateCatAtEachLoc = output_ratecatateachloc(segRateDictEst2, ratesListEst2, lenSegsEst2, dataLoc, suffix='_geopip')

output_rate_with_msa(multiAlign, rateCatAtEachLoc, dataLoc + '/msa_rate_geopip.txt')

###############################################################################
# Post-processing.

# All trees considered.
treeNames = ['true', 'phyml_1rate', 'ctmc_1rate', 'pip', 'geopip_true_start', 'geopip']
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
print 'estimated indel rates (true init)):', ratesListEst1
print 'estimated indel rates (random init)):', ratesListEst2
print 'branch length:', treeTrue.length()
print 'stationary distribution:', piProb
print 'Q:', qMat
print 'time order: phyml, ctmc, pip, geopip_true, geopip_random'
print 'time:', time1, time2, time3, time4, time5

print '### sum of square errors of all pairwise distance estimates ###'
print 'CTMC, PIP, GeoPIP (true start), GeoPIP (random start)'
print sum([(bDictCTMC1[key]-bDict[key])**2 for key in bDict.keys()])
print sum([(bDictEstPIP[key]-bDict[key])**2 for key in bDict.keys()])
print sum([(bDictEst1[key]-bDict[key])**2 for key in bDict.keys()])
print sum([(bDictEst2[key]-bDict[key])**2 for key in bDict.keys()])

print '### Weighted Robinson Foulds (wRF) of unscaled trees and true tree ###'
print json.dumps(treeDistNonScaled['true'], sort_keys=True, indent=4)

print '### Weighted Robinson Foulds (wRF) of scaled trees and true tree ###'
print json.dumps(treeDistScaled['true'], sort_keys=True, indent=4)

print '### Robinson Foulds (RF) of scaled trees trees and true tree ###'
print json.dumps(treeDistSym['true'], sort_keys=True, indent=4)
