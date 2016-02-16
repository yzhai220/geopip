# Transform between pairwise distances and trees.
# Use R code to reconstruct tree using neighor joining based on
# pairwise distances.

import numpy as np
import codecs
import dendropy
import os
from copy import deepcopy


def distmat_from_bdict(bDict):
    """
    generate a list of seqNames and a distance matrix from bDict
    input:
        bDict: dict, pair -> pairwise distance
    output:
        seqName: list, of sequence names
        dMat: matrix, of pairwise distances
    """
    nPairs = len(bDict)
    nLan = int(np.sqrt(nPairs * 2).round() + 1)
    pairsList = bDict.keys()
    set1 = set([pair[0] for pair in pairsList])
    set2 = set([pair[1] for pair in pairsList])
    seqName = list(set1.union(set2))
    seqName.sort()
    dMat = np.zeros([nLan, nLan])
    for lanPair in pairsList:
        lan1Index = seqName.index(lanPair[0])
        lan2Index = seqName.index(lanPair[1])
        dMat[lan1Index, lan2Index] = bDict[lanPair]
    dMat = dMat + dMat.T
    return seqName, dMat


def get_out_name_dist_tree_files(dataLoc, suffix=''):
    """
    generate outTreeFile location
    """
    outNameLoc = dataLoc + '/name.txt'
    outDistLoc = dataLoc + '/dist' + suffix + '.txt'
    outTreeLoc = dataLoc + '/tree' + suffix + '.txt'
    return outNameLoc, outDistLoc, outTreeLoc


def write_name(seqName, outNameLoc):
    """
    write a file containing sequence names
    """
    outNameFile = codecs.open(outNameLoc, 'w+', 'utf-8')
    strToWrite = '\t'.join(seqName)
    outNameFile.write(strToWrite)
    outNameFile.close()
    print 'Taxa names saved at %s' % (outNameLoc, )


def write_dist(dMat, outDistLoc):
    """
    write a file containing distance matrix
    """
    outFile = codecs.open(outDistLoc, 'w+', 'utf-8')
    np.savetxt(outFile, dMat)
    outFile.close()
    print 'Distance matrix saved at %s' % (outDistLoc, )


def write_dist_from_bdict(bDict, dataLoc, suffix=''):
    """
    write a file containing sequence names and a file containing distance matrix from bDict
    """
    seqName, dMat = distmat_from_bdict(bDict)
    outNameLoc, outDistLoc, outTreeLoc = get_out_name_dist_tree_files(dataLoc, suffix)
    write_name(seqName, outNameLoc)
    write_dist(dMat, outDistLoc)


def get_rscript(outNameLoc, outDistLoc, outTreeLoc, rFile='/users/vincent/googledrive/geopip/code/nj_use_r.R'):
    """
    get Rstript for constructing NJ tree from distance matrix using R code
    output:
        rCodeNj: string of code to be called by os.system()
    """
    rCodeNj = ' '.join(['Rscript', rFile, outNameLoc, outDistLoc, outTreeLoc])
    return rCodeNj


def nj_tree_from_bdict_using_r(rCodeNj, bDict, dataLoc, outTreeLoc, suffix='', rooted=True):
    """
    NJ tree from bDict, using R code, rooted (by default, using mid-point rooting) or unrooted
    """
    write_dist_from_bdict(bDict, dataLoc, suffix)
    print 'running R code to reconstruct tree using NJ'
    os.system(rCodeNj)
    tree = dendropy.Tree.get_from_path(outTreeLoc, schema='newick')
    if rooted:
        tree.reroot_at_midpoint()
    # need to read as a string again for some reason in dendropy library
    tree = dendropy.Tree.get_from_string(tree.as_newick_string(), schema="newick")
    # write rooted tree into the file, if rooted
    # treeFile = codecs.open(outTreeLoc, 'w+', 'utf-8')
    # treeFile.write(tree.as_string('newick'))
    # treeFile.close()
    return tree


def bdict_from_tree(tree, pairs):
    """
    get pairwise distance dictionary from tree
    input:
        tree: tree object
        pairs: list, of all pairs of sequences
    output:
        bDict: dict, pair -> distance
    """
    bDict = {}
    pdm = dendropy.treecalc.PatristicDistanceMatrix(deepcopy(tree))
    for pair in pairs:
        seq1Name, seq2Name = pair
        # print pdm(tree.find_node_with_taxon_label(seq1Name).taxon, tree.find_node_with_taxon_label(seq2Name).taxon)
        bDict[pair] = pdm(tree.find_node_with_taxon_label(seq1Name).taxon, tree.find_node_with_taxon_label(seq2Name).taxon)
    return bDict
