# Input and output functions for using INDELible.
#!/usr/bin/python

from io_phylip import read_phylip_msa_true


def concatenate_msa(msaFiles, removeGap=False):
    """
    Concatenate more than one MSA files.
    Args:
        msaFiles: MSA files to be concatenated, in PHYLIP format.
    Returns:
        multiAlignCon: concatenated MSA.
    """
    multiAlignCon = {}
    for msaFile in msaFiles:
        multiAlign = read_phylip_msa_true(msaFile)
        if multiAlignCon == {}:
            multiAlignCon = multiAlign
        else:
            for k, v in multiAlign.iteritems():
                vOld = multiAlignCon[k]
                vNew = vOld + v
                multiAlignCon[k] = vNew
    if removeGap:
        for k, v in multiAlignCon.iteritems():
            multiAlignCon[k] = v.replace('-', '')
    return multiAlignCon
