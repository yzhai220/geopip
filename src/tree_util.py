# Utility functions for transforming pairwise distances to and from trees.
# Use R to reconstruct tree using neighor joining based on pairwise distances.

import numpy as np
import random as rd
import codecs
import dendropy
from dendropy.calculate import treemeasure
import os


def sim_tree_newick(nLeaves, prefix='seq', bLenDist='uniform', bLenParam=[0, 0.1], bLenDecimal=4):
    """
    generate a tree with nLeaves in newick form string
    prefix: a string for each of leaves to start with, 'seq' (default)
    bLenDist: branch length distribution, 'uniform' (default), 'normal', 'beta'
    bLenParam: parames of the branch length distribution, [0, 0.1] (default)
    bLenDecimal: decimals used in the trees, 4 (default)
    """
    leavesName = [prefix+str(index) for index in range(nLeaves)]
    if bLenDist == 'uniform':
        bLens = np.random.uniform(bLenParam[0], bLenParam[1], size=2*nLeaves-2)
    elif bLenDist == 'normal':
        bLens = np.random.normal(bLenParam[0], bLenParam[1], size=2*nLeaves-2)
    elif bLenDist == 'beta':
        bLens = np.random.beta(bLenParam[0], bLenParam[1], size=2*nLeaves-2)
    else:
        print 'branch length distribution not supported'
        print 'currently only support "uniform", "normal", "beta"'
    # all branch lengths has to be positive, forcing all negative ones to 0
    bLens = [max(bLen, 0) for bLen in bLens]
    # round bLens to make it to certain decimals
    bLens = [round(bLen, bLenDecimal) for bLen in bLens]
    nodes = leavesName
    nNodes = len(nodes)
    bLenIndex = 0
    while nNodes > 1:
        nodesToCombine = rd.sample(nodes, 2)
        node1 = nodesToCombine[0]
        node2 = nodesToCombine[1]
        nodeNew = '(' + node1 + ':' + str(bLens[bLenIndex]) + ',' + node2 + ':' + str(bLens[bLenIndex+1]) + ')'
        bLenIndex += 2
        nodes.remove(node1)
        nodes.remove(node2)
        nodes.append(nodeNew)
        nNodes -= 1
    res = nodes[0] + ';'
    return res


def sim_tree_newick_fixed_edge_length(log2Leaves, bLen=0.05, prefix='seq', bLenDecimal=4):
    """
    generate a tree with 2^log2Leaves in newick form string
    prefix: a string for each of leaves to start with, 'seq' (default)
    bLenParam: parames of the branch length distribution, [0, 0.1] (default)
    bLenDecimal: decimals used in the trees, 4 (default)
    """
    nLeaves = 2**log2Leaves
    leavesName = [prefix+str(index+1) for index in range(nLeaves)]
    bLens = np.ones(2*nLeaves-2) * bLen
    nodes = leavesName
    nNodes = len(nodes)
    bLenIndex = 0
    while nNodes > 1:
        nodesToCombine = nodes[:2]
        node1 = nodesToCombine[0]
        node2 = nodesToCombine[1]
        nodeNew = '(' + node1 + ':' + str(bLens[bLenIndex]) + ',' + node2 + ':' + str(bLens[bLenIndex+1]) + ')'
        bLenIndex += 2
        nodes.remove(node1)
        nodes.remove(node2)
        nodes.append(nodeNew)
        nNodes -= 1
    res = nodes[0] + ';'
    return res


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


def get_rscript(outNameLoc, outDistLoc, outTreeLoc, rFileLoc):
    """
    get Rstript for constructing NJ tree from distance matrix using R code
    output:
        rCodeNj: string of code to be called by os.system()
    """
    rCodeNj = ' '.join(['Rscript', rFileLoc+'/nj_use_r.R', outNameLoc, outDistLoc, outTreeLoc])
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
    # Need to read as a string again for some reason in dendropy library.
    # tree = dendropy.Tree.get_from_string(tree._as_newick_string(), schema="newick")
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
    pdm = treemeasure.PatristicDistanceMatrix(tree)
    for pair in pairs:
        seq1Name, seq2Name = pair
        # print pdm(tree.find_node_with_taxon_label(seq1Name).taxon, tree.find_node_with_taxon_label(seq2Name).taxon)
        bDict[pair] = pdm(tree.find_node_with_taxon_label(seq1Name).taxon, tree.find_node_with_taxon_label(seq2Name).taxon)
    return bDict
