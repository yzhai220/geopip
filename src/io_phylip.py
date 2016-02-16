# Input and output functions for PHYLIP format.


# This function may overlap with read_phylip_infile
def read_phylip_msa_true(fileLoc):
    """
    read multiple string alignment from a file (PHYLIP format, .phylip.txt)
    input:
        fileLoc: location of the file, a string
    output:
        multiAlign: dict, seqName -> seqValue
    """
    multiAlign = {}
    f = open(fileLoc)
    # initialization for the first row
    row = f.readline()
    for row in f:
        if row == '\n':
            return multiAlign
            break
        else:
            taxaName, taxaValue = row.split()
            multiAlign[taxaName] = taxaValue


# # This function seems to fail, not used for now.
# def read_phylip_infile(fileLoc):
#     """
#     read multiple string alignment from a file (PHYLIP format, .phylip.txt)
#     input:
#         fileLoc: location of the file, a string
#     output:
#         multiAlign: dict, seqName -> seqValue
#     """
#     multiAlign = {}
#     f = open(fileLoc)
#     # initialization for the first row
#     row = f.readline()
#     # in case there are empty rows
#     while row[0] != '>':
#         row = f.readline()
#     row = row.rstrip()
#     taxaName = row[1:]
#     taxaValue = ''
#     # any other rows
#     for row in f:
#         if row[0] == '>':
#             row = row.rstrip()
#             multiAlign[taxaName] = taxaValue
#             taxaName = row[1:]
#             taxaValue = ''
#         else:
#             row = row.rstrip()
#             taxaValue += row
#     multiAlign[taxaName] = taxaValue
#     return multiAlign


def dict_write_phylip(multiAlign, outputLoc, nameLengthFixed=False, fileName='msa.phylip.txt'):
    """
    output a dict into phylip input format, with name 'msa.phylip.txt'
    input:
        multiAlign: dict, seqName -> multiple string alignment
        outputLoc: location folder of the output
    """
    outputFile = outputLoc + '/' + fileName
    f = open(outputFile, 'w')
    keysAll = multiAlign.keys()
    keysAll.sort()
    nKey = len(keysAll)
    nSite = len(multiAlign.values()[0])
    f.write(str(nKey) + '\t' + str(nSite) + '\n')
    for key in keysAll:
        keyLen = len(key)
        if nameLengthFixed:
            if keyLen < 10:
                keyStr = key + ' ' * (10-keyLen)
            else:
                keyStr = key[:10]
        else:
            keyStr = key
        seqValue = multiAlign[key]
        f.write(keyStr+'\t'+seqValue+'\n')
    f.close()


def multialign_subs_only(multiAlign):
    """
    keeps only substituions sites of a MSA
    """
    multiAlignSubs = {}
    keys = multiAlign.keys()
    values = multiAlign.values()
    msaList = zip(*values)
    msaSubsList = [msa for msa in msaList if msa.count('-') == 0]
    msaSubsList = zip(*msaSubsList)
    for index in xrange(len(keys)):
        key = keys[index]
        value = msaSubsList[index]
        value = ''.join(value)
        multiAlignSubs[key] = value
    return multiAlignSubs


def dict_write_phylip_subs_only(multiAlign, outputLoc):
    """
    output a dict into phylip input format, with name 'msaSubs.phylip.txt',
        substitions only
    input:
        multiAlign: dict, seqName -> multiple string alignment
        outputLoc: location folder of the output
    """
    multiAlignSubs = multialign_subs_only(multiAlign)
    outputFile = outputLoc + '/msaSubs.phylip.txt'
    f = open(outputFile, 'w')
    keysAll = multiAlignSubs.keys()
    keysAll.sort()
    nKey = len(keysAll)
    nSite = len(multiAlignSubs.values()[0])
    f.write(str(nKey) + '\t' + str(nSite) + '\n')
    for key in keysAll:
        keyLen = len(key)
        if keyLen < 10:
            keyStr = key + ' ' * (10-keyLen)
        else:
            keyStr = key[:10]
        seqValue = multiAlignSubs[key]
        f.write(keyStr+seqValue+'\n')
    f.close()


def dict_key_to_str(oneDict):
    newDict = {}
    for k, v in oneDict.iteritems():
        newKey = k.__str__()
        newDict[newKey] = v
    return newDict


def get_partition_from_all_pairs(pairs):
    res = []
    while len(pairs) > 0:
        pair1 = pairs[0]
        pairs.remove(pair1)
        seq1, seq2 = pair1
        for pair in pairs:
            if pair.count(seq1) == 0 and pair.count(seq2) == 0:
                pair2 = pair
                pairs.remove(pair2)
                break
        res.append([pair1, pair2])
    return res
