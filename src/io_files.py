# Input and output related functions.

import numpy as np
import codecs
import os


def get_all_files_with_same_name_in_subfolders(folder, fileName):
    """
    get all files with same name in all subfolders of one given folder
    input:
        folder: the parent folder to search through, e.g., 'tree.txt'
        fileName: a string, name of files to grab
    """
    subFolders = [folder+'/'+subFolder for subFolder in os.listdir(folder) if subFolder[0] != '.']   # only keep those subfolders which are not hidden
    allFiles = [subFolder+'/'+fileName for subFolder in subFolders]
    return allFiles


def output_list(listToWrite, outFileName):
    """
    a general function to write a list into a file
    """
    outFile = open(outFileName, 'w')
    strListToWrite = [item.__str__() for item in listToWrite]
    allToWrite = '\n'.join(strListToWrite)
    outFile.write(allToWrite)
    outFile.close()


def output_list_of_tuple(listOfTuple, outFileName):
    """
    a general function to write a list of tuples into a file
    """
    outFile = open(outFileName, 'w')
    for oneTuple in listOfTuple:
        oneTupleStr = [item.__str__() for item in oneTuple]
        oneTupleStr = '\t'.join(oneTupleStr)
        outFile.write(oneTupleStr)
        outFile.write('\n')
    outFile.close()


def output_lensegs(lenSegs, outputLoc, suffix=''):
    """
    write lenSegs into a file
    """
    outFileStr = outputLoc + '/lenSegs' + suffix + '.txt'
    output_list(lenSegs, outFileStr)


def output_ratelist(rateList, outputLoc, suffix=''):
    outFileStr = outputLoc + '/rateList' + suffix + '.txt'
    output_list_of_tuple(rateList, outFileStr)


def output_qmat(qMat, outputLoc, suffix=''):
    """
    write qMat into a file
    """
    outFileStr = outputLoc + '/qMat' + suffix + '.txt'
    outFile = codecs.open(outFileStr, 'w+', 'utf-8')
    np.savetxt(outFile, qMat)
    outFile.close()


def output_llh(llh, outputLoc, suffix=''):
    """
    write llh into a file
    """
    outFileStr = outputLoc + '/llh' + suffix + '.txt'
    outFile = open(outFileStr, 'w')
    outFile.write(llh.__str__())
    outFile.close()


def output_ratecat(rateCatList, outputLoc, suffix=''):
    """
    """
    outFileStr = outputLoc + '/rateCatList' + suffix + '.txt'
    output_list(rateCatList, outFileStr)


def ratecatlist_from_segratedict_and_ratelist(segRateDict, rateList):
    """
    calculate a list of rate category index from segRateDict and rateList
    input:
        segRateDict: dict, segId -> segRate
        rateList: list, all rates
    output:
        rateCatList: list, all
    """
    segIdAll = segRateDict.keys()
    segIdAll.sort()
    rateCatList = [rateList.index(segRateDict[segId]) for segId in segIdAll]
    return rateCatList


def ratecatateachloc_from_ratecatelist_and_lensegs(rateCatList, lenSegs):
    """
    calculate rate category for each locus from rateCatList and lenSegs
    """
    rateCatAtEachLoc = []
    for index in xrange(len(lenSegs)):
        rateCatAtEachLoc += list(np.repeat(rateCatList[index], lenSegs[index]))
    return rateCatAtEachLoc


def output_ratecatateachloc(segRateDict, rateList, lenSegs, outputLoc, suffix=''):
    rateCatList = ratecatlist_from_segratedict_and_ratelist(segRateDict, rateList)
    rateCatAtEachLoc = ratecatateachloc_from_ratecatelist_and_lensegs(rateCatList, lenSegs)
    outFileStr = outputLoc + '/rateCatAtEachLoc' + suffix + '.txt'
    strToWrite = ''.join(str(item) for item in rateCatAtEachLoc)
    f = open(outFileStr, 'w')
    f.write(strToWrite)
    f.close()
    return rateCatAtEachLoc


def output_rate_with_msa(multiAlign, rateEst, outputFile):
    values = multiAlign.values()
    rateEstStr = [rateCat.__str__() for rateCat in rateEst]
    rateEstStr = ''.join(rateEstStr)
    values.append(rateEstStr)
    output_list(values, outputFile)
