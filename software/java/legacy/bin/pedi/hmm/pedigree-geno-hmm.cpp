//////////////////////////////////////////////////
//
// Reads a pedigree and runs the following on every site
//    1) Mendelian inference
//    2) Branch-and-Bound allele assignment enumeration
//    3) Entropy calcluation for all untyped individuals in the pedigree
//
//////////////////////////////////////////////////

//#define MAX_CFGS 50
//#define DEBUG 1  // comment out to remove debugging

#include <vector>
#include <algorithm>
#include <string>
#include <set>
#include <queue>
#include <map>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <time.h>


#include "pedigree.h"
#include "pedparser.h"
#include "recombparser.h"
#include "ibd_geno_alg.h"

using namespace std;





////////////////////////
//
// ./pedigree-geno-hmm simulations/blah.rrates.1 simulations/debug.ped simulations/deleteme.ped
//







//////////////////////////////
//
//
int main (int argc, char *argv[])
{
  // read arguements
  if (argc < 3)  // argc is # argvs [1..n] with prog as argv[0]
    {
      printf ("\nUSAGE: %s [recomb file] [in pedigree file] [out pedigree file]\n\n", argv[0]);
      printf ("[recomb file]  INPUT file containing a recomb rate (double) on each line\n");
      printf ("[in pedigree file]   INPUT pedigree (genotypes ignored) in linkage format:\n");
      printf ("            FamID IndID FaID MoID Sex Aff [LOCI ...]\n");
      printf ("[out pedigree file]  OUT pedigree file name with linkage format.\n");
      printf ("\n");

      exit(-1);
    }

  const char* recomb_file = argv[1];
  const char* ped_file = argv[2];
  const char* out_ped_file = argv[3];


  printf("\nPARAMS:  recomb file = %s,\n", recomb_file);
  printf("         in ped = %s,\n", ped_file);
  printf("         out ped = %s\n\n", out_ped_file);



  // initialize random number generator
  //srand48(time(NULL));
  srand48(1);




  Pedigree pedigree;

  printf("READING input pedigree...\n");
  PedigreeParser ped_parser(1,2,0); // SNP DATA ONLY
  ped_parser.readPedigree(ped_file, pedigree);
  pedigree.setFamilyPointers();
  printf("Traversal Order:\n");
  pedigree.setTraversalOrder(true);



  printf("READING input recombination rates...\n");
  RecombParser rec_parser;
  rec_parser.openRecomb(recomb_file);
  rec_parser.readRecomb(pedigree);
  rec_parser.closeRecomb();


  //printf("recrates:  ");
  //for (int i=0; i < pedigree.recombination_rates.size();  i++)
  //{
  //  printf("  %i] %f", i, pedigree.recombination_rates[i]);
  //}
  //printf("\n");



  // set up the algorithm 
  IBDAlg ibd_alg;
  ibd_alg.setVerbose(true);


  //RecEntropyAlg rec_ent_alg;
  //rec_ent_alg,setVerbose(true);

  // Forward & Backward recursions for the IBD HMM
  //printf("Forward Recursion\n");
  //pedigree.runForwardBlockFamilyAlgorithm(ibd_alg, 1);
  //printf("Backward Recursion\n");
  pedigree.runForwardBackwardBlockFamilyAlgorithm(ibd_alg, 1);

  //pedigree.runBlockFamilyAlgorithm(rec_ent_alg, 1);

  vector<Family>::iterator fam;
  for (fam = pedigree.families.begin();  fam != pedigree.families.end();  fam++)
    {
      printf("Family %i\n", fam->family_id);
      printf("  Total Geno Entropy:  %f\n", fam->geno_entropy);
    }






  //pedigree.printInferredGenotypes(out_ped_file);  


  //printf("PRINTING Clean Linkage...\n");
  //pedigree.printSuperlinkFile(out_ped_file);


}
