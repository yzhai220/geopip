# Functions related to reading data files in fasta format.

import random as rd

from io_files import output_list


def read_fasta(fileLoc):
    """
    read multiple string alignment from a file (.fasta format)
    input:
        fileLoc: location of the file, a string
    output:
        multiAlign: dict, seqName -> seqValue
    """
    multiAlign = {}
    f = open(fileLoc)
    # initialization for the first row
    row = f.readline()
    # in case there are empty rows
    while row[0] != '>':
        row = f.readline()
    row = row.rstrip()
    taxaName = row[1:]
    taxaValue = ''
    # any other rows
    for row in f:
        if row[0] == '>':
            row = row.rstrip()
            multiAlign[taxaName] = taxaValue
            taxaName = row[1:]
            taxaValue = ''
        else:
            row = row.rstrip()
            taxaValue += row
    multiAlign[taxaName] = taxaValue
    return multiAlign


def multi_align_remove_all_empty_sites(multiAlign, step=1):
    """
    remove sites that are missing in all sequences considered
    also take a subset of remaining sites with step (by default, 1)
    """
    multiAlignNonEmptySubset = {}
    seqNameList = multiAlign.keys()
    nSeq = len(seqNameList)
    seqValueList = multiAlign.values()
    msaList = zip(*seqValueList)
    msaListNonEmpty = [msa for msa in msaList if msa.count('-') < nSeq]
    msaListNonEmpty = msaListNonEmpty[::step]
    seqValueListNew = zip(*msaListNonEmpty)
    for index in xrange(nSeq):
        seqName = seqNameList[index]
        seqValue = ''.join(seqValueListNew[index])
        multiAlignNonEmptySubset[seqName] = seqValue
    return multiAlignNonEmptySubset


def multi_align_remove_unknown_sites(multiAlign, cList):
    """
    remove sites that contains unknown symbols
    """
    multiAlignNonEmptySubset = {}
    seqNameList = multiAlign.keys()
    nSeq = len(seqNameList)
    seqValueList = multiAlign.values()
    msaList = zip(*seqValueList)
    msaListNonEmpty = [msa for msa in msaList if (msa.count('A')+msa.count('C')+msa.count('G')+msa.count('T')+msa.count('-')) == nSeq]
    seqValueListNew = zip(*msaListNonEmpty)
    for index in xrange(nSeq):
        seqName = seqNameList[index]
        seqValue = ''.join(seqValueListNew[index])
        multiAlignNonEmptySubset[seqName] = seqValue
    return multiAlignNonEmptySubset


def multi_align_remove_uncertain_sites(multiAlign, step=1, symbol='?'):
    """
    remove sites that are missing in all sequences considered
    also take a subset of remaining sites with step (by default, 1)
    """
    multiAlignNonEmptySubset = {}
    seqNameList = multiAlign.keys()
    nSeq = len(seqNameList)
    seqValueList = multiAlign.values()
    msaList = zip(*seqValueList)
    msaListNonEmpty = [msa for msa in msaList if msa.count(symbol) == 0]
    msaListNonEmpty = msaListNonEmpty[::step]
    seqValueListNew = zip(*msaListNonEmpty)
    for index in xrange(nSeq):
        seqName = seqNameList[index]
        seqValue = ''.join(seqValueListNew[index])
        multiAlignNonEmptySubset[seqName] = seqValue
    return multiAlignNonEmptySubset


def multi_align_revise_keys(multiAlign):
    """
    remove key names so that there is no space ' ', changing to ''
    """
    multiAlignNew = {}
    for k, v in multiAlign.iteritems():
        kNew = k.replace(' ', '.')
        kNew = kNew.replace(':', '.')
        kNew = kNew.replace('_', '.')
        multiAlignNew[kNew] = v
    return multiAlignNew


def multi_align_upper_case_only(multiAlign):
    """
    remove all sequences that contains elements other than cList
    """
    multiAlignUpper = multiAlign
    for k, v in multiAlign.iteritems():
        v = v.upper()
        v = v.replace('.', '-')
        multiAlignUpper[k] = v
    return multiAlignUpper


def multi_align_rna_to_dna(multiAlign):
    """
    remove all sequences that contains elements other than cList
    """
    multiAlignDna = multiAlign
    for k, v in multiAlign.iteritems():
        v = v.replace('U', 'T')
        v = v.replace('u', 't')
        multiAlignDna[k] = v
    return multiAlignDna


def multi_align_keep_only_clist(multiAlign, cList):
    """
    remove all sequences that contains elements other than cList
    """
    multiAlignClistOnly = {}
    cListExt = cList + ['-'] + ['?']
    cListExtSet = set(cListExt)
    keys = multiAlign.keys()
    for key in keys:
        vSet = set(multiAlign[key])
        if vSet.issubset(cListExtSet):
            multiAlignClistOnly[key] = multiAlign[key]
    return multiAlignClistOnly


def dict_subset(oneDict, subsetSize=10000):
    """
    get a random subset of a dict, with given size
    """
    newDict = {}
    totalSize = len(oneDict)
    if subsetSize > totalSize:
        subsetSize = totalSize
    rdIndex = rd.sample(range(totalSize), subsetSize)
    allKeys = oneDict.keys()
    subsetKeys = [allKeys[index] for index in rdIndex]
    for key in subsetKeys:
        newDict[key] = oneDict[key]
    return newDict


def dict_write_msa_no_name(multiAlign, outputLoc):
    """
    output multiAlign keys and values for visually inspection, without sequence names
    the name of the files are 'names.txt' for names of sequences, and 'msa.txt' for multple alignments
    input:
        multiAlign: dict, seqName -> multiple string alignment
        outputLoc: location folder of the output
    """
    seqNames = multiAlign.keys()
    values = multiAlign.values()
    outSeqNameFile = outputLoc + '/names.txt'
    outMsaFile = outputLoc + '/msa.txt'
    output_list(seqNames, outSeqNameFile)
    output_list(values, outMsaFile)


def dict_write_fasta(multiAlign, outputLoc, fileName='msa.fasta'):
    """
    output a dict into fasta format, with name 'msa.fasta'
    input:
        multiAlign: dict, seqName -> multiple string alignment
        outputLoc: location folder of the output
    """
    outputFile = outputLoc + '/' + fileName
    f = open(outputFile, 'w')
    keysAll = multiAlign.keys()
    keysAll.sort()
    for key in keysAll:
        seqNameStr = '>' + key + '\n'
        seqValue = multiAlign[key]
        f.write(seqNameStr)
        f.write(seqValue + '\n')
    f.close()


def multi_align_subset(multiAlign, percent=1., fixedRelativeLocation=True):
    """
    get a random subset from multiAlign, with given percentage
    input:
        multiAlign: dict, a multiple sequence alignemnt
        percent: 1 (by default), the percent of sites to be sampled
        fixedRelativeLocation: True (by default), the relative location of
            sampled sites are fixed with respect to their original location
            if False, then subsamples are shuffled
    output:
        multiAlignSub: a random sample of sites in multiAlign
    """
    multiAlignSub = {}
    seqNameList = multiAlign.keys()
    nSeq = len(seqNameList)
    seqValueList = multiAlign.values()
    msaList = zip(*seqValueList)
    nSites = len(msaList)
    nSitesSub = int(round(percent * nSites))
    if fixedRelativeLocation:
        # this sub sampling will not shuffle subsamples around
        subIndex = rd.sample(xrange(nSites), nSitesSub)
        subIndex.sort()
        msaListSub = [msaList[index] for index in subIndex]
    else:
        # note this sub sampling will shuffle subsamples around
        msaListSub = rd.sample(msaList, nSitesSub)
    seqValueListNew = zip(*msaListSub)
    for index in xrange(nSeq):
        seqName = seqNameList[index]
        seqValue = ''.join(seqValueListNew[index])
        multiAlignSub[seqName] = seqValue
    return multiAlignSub


def multi_align_subset_segs(multiAlign, nSegs=1):
    """
    get n subset from multiAlign, fixed positions
    input:
        multiAlign: dict, a multiple sequence alignemnt
        nSegs = 1
    output:
        multiAlignSegs: a list of subsets of multiAlign
    """
    multiAlignSegs = []
    seqNameList = multiAlign.keys()
    nSeq = len(seqNameList)
    seqValueList = multiAlign.values()
    msaList = zip(*seqValueList)
    nSites = len(msaList)
    nSitesSub = nSites / nSegs
    for segIndex in xrange(nSegs):
        sitesIndexStart = segIndex * nSitesSub
        sitesIndexEnd = segIndex * nSitesSub + nSitesSub
        msaListNew = msaList[sitesIndexStart:sitesIndexEnd]
        seqValueListNew = zip(*msaListNew)
        multiAlignNewSeg = {}
        for index in xrange(nSeq):
            seqName = seqNameList[index]
            seqValue = ''.join(seqValueListNew[index])
            multiAlignNewSeg[seqName] = seqValue
        multiAlignSegs.append(multiAlignNewSeg)
    return multiAlignSegs


def multi_align_subset_given_names(multiAlign, seqNames):
    res = {}
    for seqName in seqNames:
        res[seqName] = multiAlign[seqName]
    return res


def multi_align_change_key(multiAlign, keyDict):
    res = {}
    seqNames = multiAlign.keys()
    for seqName in seqNames:
        seqNameNew = keyDict[seqName]
        res[seqNameNew] = multiAlign[seqName]
    return res
