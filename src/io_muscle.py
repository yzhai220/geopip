# Input and output functions for using MUSCLE.

#!/usr/bin/python

import os


def muscle_run(dataLoc, inFName='msa_unaligned.txt.', outFName='msa.txt'):
    """
    Use muscle to get multicle sequence alignments.
    Args:
        dataLoc: data folder.
        inFName: file name of the input.
        outFName: file name of the output.
    """
    inFile = dataLoc + '/' + inFName
    outFile = dataLoc + '/' + outFName
    muscleStr = 'muscle -in ' + inFile + ' -out ' + outFile
    os.system(muscleStr)
