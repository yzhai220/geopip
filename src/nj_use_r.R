#!/usr/bin/Rscript

# Use R to contruct a NJ tree.

# Loads the command line arguments to character vector
args = commandArgs(TRUE)

# Import data files, in the order of taxa name file, distance file, and output
# tree file.
nameFile = args[1]
distFile = args[2]
outTreeFile = args[3]

# Need ape package.
ape_require = require(ape)
if (ape_require == FALSE)
{
    install.packages('ape', repos = 'http://cran.stat.sfu.ca/')
    library(ape)
}

# Load taxa name file.
name = scan(nameFile, what = 'character', sep = '\t', quiet = TRUE)

# Load distance file.
dist = read.table(distFile)
dist = as.matrix(dist)
dimnames(dist)[[1]] = name

# Contruct a tree using neighbor joining.
tree = ladderize(nj(dist))
# Set negative branch length estimates to a very small value.
tree$edge.length[tree$edge.length <= 0] = 1.e-10

# Write tree in newick form.
write.tree(tree, outTreeFile)

# Quite R.
quit()
