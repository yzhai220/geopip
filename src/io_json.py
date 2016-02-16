# Organizing input and output using json format for communication with java.

import json
import codecs
import itertools
import numpy as np
import os
import glob

from align_util import *
from pip_util import pi_from_qmat


def make_global_location_files(wd, subFolder='runs'):
    """
    make global location files:
    modelDirectory, eStepFile, parametersPath, inputLoc, outputLoc
    """
    execsLoc = wd + '/' + subFolder
    if not os.path.exists(execsLoc):
        os.makedirs(execsLoc)
    execsAll = os.listdir(execsLoc)
    execsAll = [execsOne for execsOne in execsAll if execsOne[0] != '.']
    # print execsAll
    dataLoc = execsLoc + '/' + str(len(execsAll)+1)
    inputLoc = dataLoc + '/input'
    outputLoc = dataLoc + '/output'
    if not os.path.exists(inputLoc):
        os.makedirs(inputLoc)
    if not os.path.exists(outputLoc):
        os.makedirs(outputLoc)
    eStepFile = dataLoc + '/Estepfile.json'
    parametersPath = dataLoc + '/input.json'
    return eStepFile, parametersPath, inputLoc, outputLoc, dataLoc


def create_states_json(cList, fLoc):
    """
    create a 'states.json' file used by java code as the model
    to the fLoc folder
    """
    fName = fLoc + '/states.json'
    outFile = codecs.open(fName, 'w', 'utf-8')
    json.dump(cList, outFile)
    outFile.close()
    return


def create_features_json(cList, fLoc):
    """
    create a 'feature.json' file used by java code as the model
    to the fLoc folder
    this one assumes GTR model
    """
    res = []
    fName = fLoc + '/features.json'
    outFile = codecs.open(fName, 'w', 'utf-8')
    states = cList
    nCh = len(cList)
    bases = ['base=' + ch for ch in cList]
    pairsIter = itertools.combinations(cList, 2)
    pairs = [list(pair) for pair in pairsIter]
    nPair = len(pairs)
    transitions = ['isTransition=' + ''.join(pair) for pair in pairs]
    res += [{'features': [bases[i]], 'states':[states[i]]} for i in range(nCh)]
    res += [{'features': [transitions[i]], 'states':pairs[i]} for i in range(nPair)]
    json.dump(res, outFile, indent=4)
    outFile.close()
    return


def write_align_2lan(lanName1, lanName2, t, alignPath, alignRes):
    """
    create two txt files:
    lanName1_lanName2.tree.txt:
    tree file, with t splited into half for each, (lanName1:t/2, lanName2:t/2);
    lanName1_lanName2.align.txt:
    alignment file, with format
    >lanName1
    ACGTAAAA....
    >lanName2
    AACTACGC....
    """
    names = [lanName1, lanName2]
    # names = sorted(names)
    nameLinked = '_'.join(names)
    if not os.path.exists(alignPath):
        os.makedirs(alignPath)
    treeFilePath = alignPath + '/' + nameLinked + '.tree.txt'
    alignFilePath = alignPath + '/' + nameLinked + '.align.txt'
    # alignment file writing
    alignF = open(alignFilePath, 'w')
    alignPre1 = '>' + lanName1 + '\n'
    alignPre2 = '>' + lanName2 + '\n'
    alignReduced = align_reduce(alignRes)
    str1, str2 = align_str(alignReduced)
    alignF.write(alignPre1)
    alignF.write(str1)
    alignF.write('\n')
    alignF.write(alignPre2)
    alignF.write(str2)
    alignF.write('\n')
    alignF.close()
    # tree file writing
    treeF = open(treeFilePath, 'w')
    # NOTE: THERE MIGHT BE A BETTER WAY TO INPUT BOTH STRING AND NUMBER
    tree = '(' + lanName1 + ':' + str(t/2) + ',' + lanName2 + ':' + str(t/2) + ');'
    treeF.write(tree)
    treeF.close()
    return


def write_align_2lan_tree_only(lanName1, lanName2, t, alignPath):
    """
    create one tree txt files:
    lanName1_lanName2.tree.txt:
    tree file, with t splited into half for each, (lanName1:t/2, lanName2:t/2);
    """
    names = [lanName1, lanName2]
    # names = sorted(names)
    nameLinked = '_'.join(names)
    if not os.path.exists(alignPath):
        os.makedirs(alignPath)
    treeFilePath = alignPath + '/' + nameLinked + '.tree.txt'
    # tree file writing
    treeF = open(treeFilePath, 'w')
    # NOTE: THERE MIGHT BE A BETTER WAY TO INPUT BOTH STRING AND NUMBER
    tree = '(' + lanName1 + ':' + str(t/2) + ',' + lanName2 + ':' + str(t/2) + ');'
    treeF.write(tree)
    treeF.close()
    return


def write_align(pairs, aligns, ts, alignPath):
    """
    create all txt files in pairs:
    lanName1_lanName2.tree.txt:
    tree file, with t splited into half for each, (lanName1:t/2, lanName2:t/2);
    lanName1_lanName2.align.txt:
    alignment file, with format
    >lanName1
    ACGTAAAA....
    >lanName2
    AACTACGC....
    """
    for pair in pairs:
        lanName1 = pair[0]
        lanName2 = pair[1]
        alignRes = aligns[(lanName1, lanName2)]
        t = ts[(lanName1, lanName2)]
        write_align_2lan(lanName1, lanName2, t, alignPath, alignRes)
    return


# TODO: may combine this one with previous one later
def write_align_tree_only(pairs, ts, alignPath):
    """
    create only tree txt files in pairs when alignment may stay the same:
    lanName1_lanName2.tree.txt:
    tree file, with t splited into half for each, (lanName1:t/2, lanName2:t/2);
    """
    for pair in pairs:
        lanName1 = pair[0]
        lanName2 = pair[1]
        # alignRes = aligns[(lanName1, lanName2)]
        t = ts[(lanName1, lanName2)]
        write_align_2lan_tree_only(lanName1, lanName2, t, alignPath)
    return


def dict_write_align_fasta(multiAlign, outputFile):
    """
    output a dict into fasta format, with name fileName
    input:
        multiAlign: dict, seqName -> multiple string alignment
        outputLoc: location folder of the output
    """
    f = open(outputFile, 'w')
    keysAll = multiAlign.keys()
    keysAll.sort()
    for key in keysAll:
        seqNameStr = '>' + key + '\n'
        seqValue = multiAlign[key]
        f.write(seqNameStr)
        f.write(seqValue + '\n')
    f.close()


def write_tree(tree, outTreeFile):
    f = open(outTreeFile, 'w')
    treeStr = tree.as_string('newick')
    treeStrStartIndex = treeStr.index('(')
    treeStr = treeStr[treeStrStartIndex:]
    f.write(treeStr)
    f.close()


def get_java_jars(javaDirectory):
    """
    Get all jar files needed for running Java code to estimate rate matrix Q.
    """
    # Get directories.
    dir_jars = javaDirectory + '/jars'
    jars = glob.glob(dir_jars + '/*.jar')
    legacy_str = javaDirectory + '/legacy/bin'
    jars_str = legacy_str + ':' + ':'.join(jars)
    javaJars = 'java -Xmx8g -cp  ' + jars_str
    return javaJars


def java_mstep(modelDirectory, eStepFile, javaDirectory):
    """
    generate command line code for java M step
    """
    javaJars = get_java_jars(javaDirectory)
    res = javaJars
    res += ' conifer.ml.main.MMain'
    res += ' -modelDirectory ' + modelDirectory    # input: model directory
    res += ' -EStepfiles ' + eStepFile    # input: E step files
    # res += ' -execDir ' + outputDirectory
    return res


def java_estep(treeFile, alignmentFile, parametersPath, branchOutput, rootOutput, likelihoodFile, javaDirectory):
    """
    generate command line code for java E step
    """
    javaJars = get_java_jars(javaDirectory)
    res = javaJars
    res += ' conifer.ml.ComputePosteriorMain'
    res += ' -treeFile ' + treeFile    # input: tree file, txt
    res += ' -alignmentFile ' + alignmentFile    # input: alignment file, txt
    res += ' -parametersPath ' + parametersPath    # input: params file, json
    res += ' -branchOutput ' + branchOutput    # output: branchs, into a folder
    res += ' -rootOutput ' + rootOutput    # output: root, json
    res += ' -likelihood ' + likelihoodFile    # output: likelihood, txt
    res += ' -outputCTMCStatistics true'
    return res


def param_for_estep(qMat, rate=1.0, categoryPriors=1.0):
    """
    generate the input.json file required for the E step
    contains rate matrix and categories of rates
    NOTE: THIS FUNCTION CURRENTLY ONLY WORK FOR ONE RATE
          DIFFERENT FROM R CODE INPUT
    """
    res = {}
    # nCharType = qMat.shape[0]
    # nRate = length(rate)
    # nCate = length(categoryPriors)
    piProb = pi_from_qmat(qMat)
    res['categoryPriors'] = [1]
    res['observationErrorProbability'] = 0
    res['stationaryDistributions'] = [piProb]
    res['rateMatrices'] = [qMat]
    return res


class NumPyArangeEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()    # or map(int, obj)
        return json.JSONEncoder.default(self, obj)


def write_param_for_estep(qMat, parametersPath, rate=1.0, categoryPriors=1.0):
    """
    write the parameter file for Estep, i.e., input.json in R code
    """
    res = param_for_estep(qMat, rate=1.0, categoryPriors=1.0)
    outFile = codecs.open(parametersPath, 'w', 'utf-8')
    json.dump(res, outFile, cls=NumPyArangeEncoder)
    outFile.close()
    return


def summarize_estep_one_pair(listFiles, rootFile):
    """
    summarize results from E step files into one file to be used for M step
    including holding times, wait times, and root distribution
    """
    fRoot = open(rootFile)
    nRoot = np.array(json.load(fRoot)[0])
    nCharType = len(nRoot)
    tranMat = np.zeros((nCharType, nCharType))
    for listFile in listFiles:
        fTem = open(listFile)
        tranMat += np.array(json.load(fTem)[0])
    hTime = tranMat.diagonal()
    return tranMat, hTime, nRoot


def summarize_estep_all_pairs(listFiles, rootFiles, nCharType):
    """
    summarize results from E step files of all pairs
    """
    nPair = len(rootFiles)
    hTime = np.zeros(nCharType)
    nRoot = np.zeros(nCharType)
    tranMat = np.zeros((nCharType, nCharType))
    for pairIndex in range(nPair):
        tranMatTem, hTimeTem, nRootTem = summarize_estep_one_pair(listFiles[pairIndex], rootFiles[pairIndex])
        tranMat += tranMatTem
        hTime += hTimeTem
        nRoot += nRootTem
    return tranMat, hTime, nRoot


def write_estepfiles_for_mstep(listFiles, rootFiles, eStepFile, nCharType, cList):
    """
    write a EStepFile.json for M step
    """
    tranMat, hTime, nRoot = summarize_estep_all_pairs(listFiles, rootFiles, nCharType)
    res = []
    for ch in cList:
        chIndex = cList.index(ch)
        res += [{'states': [ch], 'type':'INIT', 'value':nRoot[chIndex]}]
        res += [{'states': [ch], 'type':'HOLD', 'value':hTime[chIndex]}]
    pairsIter = itertools.permutations(cList, 2)
    pairs = [list(pair) for pair in pairsIter]
    for pair in pairs:
        ch1 = pair[0]
        ch2 = pair[1]
        ch1Index = cList.index(ch1)
        ch2Index = cList.index(ch2)
        res += [{'states': [ch1, ch2], 'type':'TRANS', 'value':tranMat[ch1Index, ch2Index]}]
    # print res
    outFile = codecs.open(eStepFile, 'w+', 'utf-8')
    json.dump(res, outFile, indent=4)
    outFile.close()
    return


def find_newest_exec(execsLoc):
    """
    find the latest exec in the execs folder
    return the lastest folder
    """
    execsAll = os.listdir(execsLoc)
    execsAllNum = [int(execName.split('.')[0]) for execName in execsAll if execName.endswith('exec')]
    execsAllNum.sort()
    execNewNum = execsAllNum[-1]
    execNew = str(execNewNum) + '.exec'
    execNewFullPath = execsLoc + '/' + execNew
    return execNewFullPath


def find_newest_exec_plus_one(execsLoc):
    """
    find the latest exec in the execs folder
    return the lastest folder
    """
    execsAll = os.listdir(execsLoc)
    execsAllNum = [int(execName.split('.')[0]) for execName in execsAll if execName.endswith('exec')]
    execsAllNum.sort()
    execNewNum = execsAllNum[-1]
    execNewNum += 1
    execNew = str(execNewNum) + '.exec'
    execNewFullPath = execsLoc + '/' + execNew
    return execNewFullPath


def summarize_mstep(execNewFullPath):
    """
    summarize M Step results
    """
    qMatPath = execNewFullPath + '/learned-rate-matrix.json'
    qMatFile = open(qMatPath)
    qMat = np.array(json.load(qMatFile))
    piProbPath = execNewFullPath + '/learned-stat-dist.json'
    piProbFile = open(piProbPath)
    piProb = np.array(json.load(piProbFile))
    wPath = execNewFullPath + '/learned-weights.json'
    wFile = open(wPath)
    w = json.load(wFile)
    return qMat, piProb, w


def write_param_estep_from_mstep(execNewFullPath, parametersPath, rate=1.0, categoryPriors=1.0):
    """
    find the Q estimate from M Step results and move it to the E Step input
    folder to be used in E Step
    """
    qMat, piProb, w = summarize_mstep(execNewFullPath)
    write_param_for_estep(qMat, parametersPath, rate=1.0, categoryPriors=1.0)
    return qMat, piProb


def get_treefiles_alignfiles(inputLoc):
    """
    get all tree files and align files from folder inputLoc
    """
    fAll = os.listdir(inputLoc)
    trees = [f for f in fAll if f.endswith('tree.txt')]
    trees.sort()
    aligns = [f for f in fAll if f.endswith('align.txt')]
    aligns.sort()
    lanPairNames = [tree.split('.')[0] for tree in trees]
    treeFiles = [inputLoc + '/' + tree for tree in trees]
    alignFiles = [inputLoc + '/' + align for align in aligns]
    treeAlignFiles = zip(treeFiles, alignFiles)
    treeAlignDict = dict(zip(lanPairNames, treeAlignFiles))
    return treeAlignDict


def java_estep_batch(treeAlignDict, outputLoc, parametersPath, javaDirectory):
    """
    return a batch of java E step codes for each pair of languages
    """
    res = []
    for lanPairName, treeAlignFile in treeAlignDict.iteritems():
        treeFile = treeAlignFile[0]
        alignmentFile = treeAlignFile[1]
        outputFolderName = outputLoc + '/' + lanPairName
        if not os.path.exists(outputFolderName):
            os.makedirs(outputFolderName)
        branchOutput = outputFolderName
        rootOutput = outputFolderName + '/' + 'root.txt'
        likelihoodFile = outputFolderName + '/llh.txt'
        javacodeE = java_estep(treeFile, alignmentFile, parametersPath, branchOutput, rootOutput, likelihoodFile, javaDirectory)
        res.append(javacodeE)
    return res


def filepath_from_output_folder(outputLoc):
    """
    get all listFiles (ctmc statistics), rootFiles (root stationary est), and
    logLililihood file names from output folder
    """
    subFolders = os.listdir(outputLoc)
    subFolders = [subFolder for subFolder in subFolders if not subFolder.startswith('.')]
    listFiles = []
    rootFiles = []
    llhFiles = []
    for subFolder in subFolders:
        subFolderPath = outputLoc + '/' + subFolder
        allFiles = os.listdir(subFolderPath)
        listFilesNew = [subFolderPath + '/' + oneFile for oneFile in allFiles if oneFile.endswith('ctmcStatistics.json')]
        rootFileNew = subFolderPath + '/root.txt'
        llhFileNew = subFolderPath + '/llh.txt'
        listFiles.append(listFilesNew)
        rootFiles.append(rootFileNew)
        llhFiles.append(llhFileNew)
    return listFiles, rootFiles, llhFiles


def llh_from_llhfiles(llhFiles):
    """
    calculate loglikelihood from all loglikelihood files
    """
    llh = 0
    for llhFile in llhFiles:
        fileTem = open(llhFile)
        llhNew = float(fileTem.read())
        llh += llhNew
    return llh


def write_multi_alignment(multiAlignSeqs, dataLoc):
    """
    write multiple alignment into file
    '*' is the symble used to separate segments
    '>' is the symble used to separate taxa name with taxa value
    """
    fileName = dataLoc + '/multiAlign.txt'
    f = codecs.open(fileName, 'w+', 'utf-8')
    for ka, va in multiAlignSeqs.iteritems():
        f.write(ka)
        f.write('>')
        f.write(va)
        f.write('\n')
    f.close()
    print 'Multiple alignment saved at %s' % (fileName, )
    return
