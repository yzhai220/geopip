#
# Makefile params
#
 
 
all: pedigree-bnb pedigree-entropy


#OPT = -O3
#OPT = -O0 -pg
OPT = -O0 -g



pedigree-bnb: SymbolTree.o pedigree-bnb.cpp pedigree.h familyalgorithm.h mendelian_inference_alg.h bnb_alg.h genoparser.h pedparser.h Makefile
	g++ $(OPT) -Wall -o pedigree-bnb SymbolTree.o pedigree-bnb.cpp

pedigree-hmm: SymbolTree.o pedigree-hmm.cpp pedigree.h familyalgorithm.h ibd_set_forest.h ibd_alg.h genoparser.h pedparser.h recombparser.h Makefile
	g++ $(OPT) -Wall -o pedigree-hmm SymbolTree.o pedigree-hmm.cpp

pedigree-geno: SymbolTree.o pedigree-geno-hmm.cpp pedigree.h familyalgorithm.h ibd_set_forest.h ibd_geno_alg.h genoparser.h pedparser.h recombparser.h Makefile
	g++ $(OPT) -Wall -o pedigree-geno-hmm SymbolTree.o pedigree-geno-hmm.cpp


SymbolTree.o: SymbolTree.cpp SymbolTree.h Makefile
	g++ $(OPT) -Wall -c SymbolTree.cpp





#
# cleaning params
#
clean:
	rm ./*.o ./tag  ./core ./*~ 
