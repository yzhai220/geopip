Supplementary Source Code and Examples for 

“A Poissonian model of indel rate variation for phylogenetic tree inference”

by Yongliang Zhai and Alexandre Bouchard-Côté.

Department of Statistics, the University of British Columbia, Vancouver, Canada.

Email: y.zhai@stat.ubc.ca, bouchard@stat.ubc.ca

################################################################################

Python, Java and R are required for running the scripts. Make sure that you have

1. Java SE Development Kit 8 (JDK) installed.

2. R installed, and Rscript callable from the command line directly. You may need to install the newest version of R before running the scripts.

3. R package 'ape' installed.

4. Python package 'dendropy' installed.

Software PhyML, INDELible and MUSCLE are included with this release for comparison purpose. See http://www.atgc-montpellier.fr/phyml/, http://abacus.gene.ucl.ac.uk/software/indelible/ and http://www.drive5.com/muscle/ for more information.

################################################################################

To run a script in the command line, for example, GeoPIP simulation studies, navigate to the /script subfolder and then type in command line:

     python geopip_test.py

The results will be generated in the /result folder. The print.txt file summarizes the results. 

################################################################################

List of subfolders:

1. script: examples of simulation studies and data analysis used in the paper.

2. src: all source codes in Python and R, developed for this paper.

3. software: all source codes in Java, developed for our previous work, and INDELible and MUSCLE.

4. result: folder for results output.

5. model: substitution model configuration files for Java code. The GTR and HKY85 model are provided.

6. data: Molluscan data used in the paper. This data set was downloaded from http://www.rna.icmb.utexas.edu/SIM/4D/Mollusk/alignment.gb, and then changed to the fasta format.
