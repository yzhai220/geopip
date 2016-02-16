
///////////////////////////////////////////
//
//  Pedigree data structures:
//      Marker
//      Individual
//      Family
//  Pedigree methods
//      read_allele
//      read_pedigree_file
//      
//      

#ifndef PEDIGREE
#define PEDIGREE



#include <vector>
#include <algorithm>
#include <string>
#include <set>
#include <map>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <time.h>


//ASSERTS
#define DEBUGASSERT

#ifndef DEBUGASSERT
#define assert(x)
#else
#define assert(x)	      \
  if (! (x))		      \
    { \
       cout << "ERROR!! Assert " << #x << " failed\n"; \
       cout << " on line " << __LINE__  << "\n"; \
       cout << " in file " << __FILE__ << "\n";  \
    }
#endif


#include "indicator.h"
#include "familyalgorithm.h"
#include "state_space.h"
#include "ibd_set_forest.h"







using namespace std;

#ifndef UNKNOWN_LITERAL
#define UNKNOWN_LITERAL -33
#endif


#define UNKNOWN 0

#define MALE 1
#define FEMALE 2
#define SAME_AS 10

#define AFFECTED 2
#define UNAFFECTED 1

#define FIRST_ALLELE 1
#define SECOND_ALLELE 2


#define INVALID -1
#define ERROR -55


#define INVALID_ALLELE -22


#define CHECK_F 7  // debugging output for a particular family
#define CHECK_I 7  // debugging output for a particular individual

#define DEBUG 1  // comment out to remove debugging




const double log_half = log(0.5);


typedef int Gender;
typedef int Affection;




class Pedigree;
class Family;
class Individual;
class Marker;

class OrderItem;
class FounderChildDLists;
class FounderPairDLists;
class DistinctAlleles;

class FamilyAlgorithm;




double mylog (double x);
int intpow (int x, int y);




void haplotype_to_int(long int& haplotype1, long int& haplotype2, const vector<Marker>& genotype_array, int start_snp, int upper_bound_snp);
long int haplotype_to_int(vector<int> allele_array);
void int_to_haplotype(vector<int>& allele_array, const long int haplotype, const int number_of_snps, const bool doPrint, FILE* fp);








class Marker
{
 public:
  Marker(){
    allele1 = UNKNOWN;
    allele2 = UNKNOWN;
  }

  // make sure alleles are sorted when storing them as markers
  Marker(int a, int b);

  Marker(const Marker& m){
    allele1 = m.allele1;
    allele2 = m.allele2;
  }
  const Marker& operator= (const Marker& m) {
    allele1 = m.allele1;
    allele2 = m.allele2;
    return *this;
  }

  // insert an allele in the case that there can be at most one copy of each type of allele
  int insertDistinctAllele(int a);

  int allele1;
  int allele2;
};


Marker::Marker(int a, int b)
{
  if (a <= b)
    {
      allele1 = a;
      allele2 = b;
      if (allele1 == UNKNOWN && allele2 == FIRST_ALLELE)
	{
	  allele1 = b;
	  allele2 = a;
	}
    }
  else
    {
      allele1 = b;
      allele2 = a;	
      if (allele1 == UNKNOWN && allele2 == FIRST_ALLELE)
	{
	  allele1 = a;
	  allele2 = b;
	}
    }
}




// insert an allele in the case that there can be at most one copy of each type of allele
int Marker::insertDistinctAllele(int a)
{
  if (a <= UNKNOWN)
    {
      return INVALID_ALLELE;
    }
  if (allele1 != UNKNOWN && allele2 != UNKNOWN 
      && (a != allele1 || a != allele2))
    {
      return INVALID_ALLELE;
    }
  if (a == allele1 || a == allele2)
    return 0;
    
  if (allele1 == UNKNOWN && allele2 == UNKNOWN)
    {
      if (a == FIRST_ALLELE)
	allele1 = a;
      if (a == SECOND_ALLELE)
	allele2 = a;
      return 1;
    }
  else
    {
      if (allele1 == UNKNOWN)
	allele1 = a;
      else
	allele2 = a;
      if (allele1 > allele2)
	return ERROR;
      return 1;
    }
  return 0;
}








class Individual
{
 public:
  Individual();
  Individual(int fid, int id, int f, int m, Gender s, Affection a, long p, int l);
  Individual(const Individual& i);
  const Individual& operator= (const Individual& i);

  ~Individual(){
    if (likelihood != NULL)
      delete [] likelihood;
    if (pair_likelihood != NULL)
      delete [] pair_likelihood;
  }

  bool operator== (const Individual& indiv) {
    if (family_id == indiv.family_id)
      return this->individual_id == indiv.individual_id;
    else
      return 0;
  }

  void setIndivIndex(int i){
    individual_index = i;
  }
  void setNumHaps(int numofhaps){
    num_haps = numofhaps;
  }


  bool haplotypeFits(long haplotype, 
		    const unsigned int start_snp, 
		    const unsigned int upper_bound_snp);


  void initPairLikelihood(int numofhaps, double init);

  void set_same_parental_source(int locus, unsigned int same_index);
  void set_value_parental_source(int locus, int new_value);
  void flatten_parental_source();
  int parental_source_root(int locus);


  int family_id;
  int individual_id;
  int individual_index;
  int father_id;
  int mother_id;
  Individual* mother;
  Individual* father;
  Gender sex;
  Affection affection;
  int isTyped;

  int order;

  long file_position;
  int line_number;

  // True if the parental source and genotypes
  // are such that any haplotype can fit in this individual.
  // Actually under the most stringent possible conditions regarding PS settings,
  // i.e. there could be other conditions where one could set any_haplotype.
  int any_haplotype;


  vector<Individual*> children; 
  vector<Marker> genotypes;
  vector<Marker> inferred_genotypes;


  // one source for each marker
  // A glorified set data structure
  // 0 -> unknown
  // 1 -> allele a (a < b) from father
  // 2 -> allele a (a < b) from mother
  // SAME_AS -> parental_source[j] = parental_source[offset]
  vector<int> parental_source; 
  int isSymmetric;
  bool hetero_unknown;

  // Haplotype vector.
  // Entries come in pairs.
  // All haplotypes must be valid [0..2^S-1] 
  // where S is the number of SNPs.
  // In the case that isSymmetric == true, the 
  // first haplotype from the pair is the one
  // from the father.
  vector<long> haplotype_pairs;

  vector<long> sampled_haplotype_pairs;
  map<long int, int> sampled_counts;
  long int sample_max_key;


  long int haplotype1;
  long int haplotype2;
  int haplotype1_source; // MALE, FEMALE, or UNKNOWN

  double log_likelihood;


  // allocated only in people who are  not founders or children of founders
  //vector<double> likelihood; // normalized probabilities
  double* likelihood; // normalized probabilities
  int num_haps;

  // allocated only in children of founders
  // an upper-right matrix (lower right is zero'ed)
  //vector<vector<double> > pair_likelihood;  // normalized probabilities
  double* pair_likelihood;  // normalized probabilities

  // the sum of values probabilities in either likelihood or pair_likelihood
  // (if likelihood or pair_likelihood contain log(p) where p is a probability
  // then sum contains the sum of p's).
  double sum;


  bool aff_hap;
};











class OrderItem
{
 public:
  OrderItem(){
    individual_id = UNKNOWN;
    individual_index = INVALID;
    individual = NULL;
    depth = INVALID;

  }
  OrderItem(int id){
    individual_id = id;
    individual_index = INVALID;
    individual = NULL;
    depth = INVALID;

  }
  OrderItem(const OrderItem& o){
    individual_id = o.individual_id;
    individual_index = o.individual_index;
    individual = o.individual;
    depth = o.depth;
  }
  const OrderItem& operator= (const OrderItem& o) {
    individual_id = o.individual_id;
    individual_index = o.individual_index;
    individual = o.individual;
    depth = o.depth;
    return *this;
  }

  int individual_id;
  int individual_index;
  Individual* individual;
  int depth;
};
struct order_lt
{
  bool operator()(const OrderItem i, const OrderItem j) const
  {
    return i.depth > j.depth;
  }
};





class FounderChildDLists
{
public:
  FounderChildDLists(Individual* c)
  {
    founder_child = c;
  }

  Individual* founder_child;
  set<int> descent_list;
};


class FounderPairDLists
{
public:
  FounderPairDLists(Individual* i, Individual* j)
  {
    founding_pair.push_back(i);
    founding_pair.push_back(j);
  }
  ~FounderPairDLists()
  {
    vector<FounderChildDLists*>::iterator cdlist;
    for (cdlist = full_children.begin();  cdlist != full_children.end();  cdlist++)
      {
	delete (*cdlist);
      }
  }

  vector<Individual*> founding_pair;  // cardinality: 2
  vector<FounderChildDLists*> full_children; // cardinality: # of full siblings

};





class Family
{
public:
  Family(int id){ 
    family_id = id; 
    isCompatible = false;
    log_likelihood = -999999999999.9;
    geno_backward = NULL;
    geno_backward_current = NULL;
    haplo_backward = NULL;
    haplo_backward_current = NULL;
    emission_prob = NULL;
    forest = NULL;
    alt_sum = 0.0;
    null_sum = 0.0;
  }
  ~Family(){
    if (forest != NULL){
      delete forest;
      forest = NULL;
    }
  }

  bool operator== (const Family& fam) {
    return this->family_id == fam.family_id;
  } 
  bool operator== (vector<Family>::iterator fam) {
    return this->family_id == fam->family_id;
  }


  void setIndices(){
    for(unsigned int i = 0;  i < members.size();  i++){
      members[i]->setIndivIndex(i);
    }
  }
  void setPedigree(Pedigree* _pedigree);
  void setFamilyPointers();
  void setTraversalOrder(bool verbose);
  void dfs(Individual* parent, int depth);

  void setNumMeioses();
  void setNumAltMeioses();


  void clearGenotypes(){
    for (unsigned int i = 0;  i < members.size(); i++)
      members[i]->genotypes.clear();
  }



  void printInconsistencies();



  void clearCandidateHaplotypes(){
    if (!candidate_haplotypes.empty())
      candidate_haplotypes.clear();
  }
  void setCandidateHaplotype(long h){
    // only insert h if it isn't already in the list
    pair<set<long>::iterator, bool> p;
    p = candidate_haplotypes.insert(h);
  }
  bool isCandidateHaplotype(long h){
    // returns true if the haplotype is found
    set<long>::iterator found;
    found = candidate_haplotypes.find(h);
    if (found != candidate_haplotypes.end())
      return true;
    return false;
  }
  void printCandidateHaplotypes(){
    printf("Candidate haplotypes for fam %i:\n      ", family_id);
    set<long>::iterator ch;
    for (ch = candidate_haplotypes.begin();  
	 ch != candidate_haplotypes.end();  
	 ch++)
      {
	printf("%li ", *ch);
      }
    printf("\n");
  }
  void copyCandidateHaplotypesInto(vector<Family>::iterator fam){
    fam->candidate_haplotypes = candidate_haplotypes;
  }
  int getCandidateHaplotypesSize(){
    return candidate_haplotypes.size();
  }

  void clearConsistencyVars(){
    null_sum = 0;
    alt_sum = 0;

    for (unsigned int i = 0;  i < geno_dynprog.size(); i++)
      delete geno_dynprog[i];
    geno_dynprog.clear();
  }


  vector<Individual*>::iterator find_individual(vector<Individual*>::iterator curr, 
						vector<Individual*>::iterator end, 
						int individual_id);





  int family_id;
  vector<Individual*> members;

  vector<Individual*> founder; // index of the child in the vector of individuals
  vector<Individual*> nonfounder; // index of the child in the vector of individuals
  vector<Individual*> typed; // index of the child in the vector of individuals
  vector<Individual*> untyped; // index of the child in the vector of individuals


  set<long> candidate_haplotypes;

  vector<DistinctAlleles>* pedigree_marker_alleles;

  vector<OrderItem> order;  // bottom up order
  vector<FounderPairDLists*> descent_lists;

  bool isCompatible;
  double log_likelihood; // transmission likelihood

  Pedigree* pedigree;


  // HMM Dynamic Programming Data Structures

  vector<GenoProb*> geno_dynprog;
  vector<HaploProb*> haplo_dynprog;


  double geno_log_prob_data;
  double haplo_log_prob_data;

  double geno_entropy;
  double haplo_entropy;

  GenoProb* geno_backward;
  GenoProb* geno_backward_current;
  HaploProb* haplo_backward;
  HaploProb* haplo_backward_current;

  GenoProb* emission_prob; // current emission on backward pass

  vector<IBDSetForest> curr_inheritance; // current inheritance

  // 'global' results for the forward-backward alg.
  //vector<vector<IBDSetForest> > inheritance;
  vector<vector<double> > marginal_indicators;

  // vector continaing weights for inheritance paths
  // indexes: site, indicator index
  vector<vector<Indicator> > indicator_weights;


  long rand_inherit;
  IBDSetForest* forest;

  double null_sum;
  double alt_sum;


  double consistency_log_likelihood;

  double log_sample_prob;

  int num_meioses;
  int num_alt_meioses;
};





class DistinctAlleles
{
 public:
  void clear(){
    if (!alleles.empty())
      alleles.clear();
  }
  void insert(int a){
    // only insert a if it isn't already in the list
    pair<set<int>::iterator, bool> p;
    p = alleles.insert(a);
  }
  bool has(int a){
    // returns true if the haplotype is found
    set<int>::iterator found;
    found = alleles.find(a);
    if (found != alleles.end())
      return true;
    return false;
  }
  void print(){
    set<int>::iterator ch;
    for (ch = alleles.begin();  
	 ch != alleles.end();  
	 ch++)
      {
	printf("%i ", *ch);
      }
    printf("\n");
  }
  int size(){
    return alleles.size();
  }

 protected:
  set<int> alleles;
};







class Pedigree
{
 public:
  Pedigree(){};
  ~Pedigree();


  void setAlleleBounds(int _smallest_allele, int _largest_allele);
  void setFamilyPointers();
  void setTraversalOrder(bool verbose);

  void setNumMeioses();
  void setNumAltMeioses();

  
  void runFamilyAlgorithm(FamilyAlgorithm& alg);
  void runForwardBlockFamilyAlgorithm(FamilyAlgorithm& alg, int block_size);
  void runBackwardBlockFamilyAlgorithm(FamilyAlgorithm& alg, int block_size);
  void runForwardBackwardBlockFamilyAlgorithm(FamilyAlgorithm& alg, int block_size);

  void printSuperlinkFile(const char* hap_file);
  void printInferredGenotypes(const char* hap_file);
  void printInferredGenotypes();
  void printGenotypes();
  void printInferredMultilocusHaps(const char* hap_file, int block_size);

  void clearGenotypes();

  void clearConsistencyVars();



  vector<Family> families;
  vector<DistinctAlleles> marker_alleles;

  vector<double> recombination_rates;


 protected:
  int smallest_allele;
  int largest_allele;


};




Pedigree::~Pedigree()
{
  // For each family
  vector<Family>::iterator fam;
  for (fam = families.begin();  fam != families.end();  fam++)
    {
      if(fam->forest != NULL)
	{
	  delete fam->forest;
	  fam->forest = NULL;
	}
      // For each person in the family
      vector<Individual*>::iterator person;
      for (person = fam->members.begin();  person != fam->members.end();  person++)
	{
	  delete (*person);
	}
	
      vector<GenoProb*>::iterator marker;
      for (marker = fam->geno_dynprog.begin();  marker != fam->geno_dynprog.end();  marker++)
	{
	  delete (*marker);
	}

      // Delete the descent_lists
      vector<FounderPairDLists*>::iterator fdlist;
      for(fdlist = fam->descent_lists.begin();  fdlist != fam->descent_lists.end();  fdlist++)
	{
	  delete (*fdlist);
	}
    }  
  families.clear();
}


void Pedigree::setAlleleBounds(int _smallest_allele, int _largest_allele)
{
  if (_smallest_allele <= UNKNOWN)
    printf("FATAL ERROR in GenoParser(%i,%i): smallest allele <= %i\n", _smallest_allele, _largest_allele, UNKNOWN);
  if (_smallest_allele > _largest_allele)
    printf("FATAL ERROR in GenoParser(%i,%i): smallest allele > largest allele\n", _smallest_allele, _largest_allele);
  
  smallest_allele = _smallest_allele;
  largest_allele = _largest_allele;
}



void Pedigree::clearGenotypes()
{
  vector<Family>::iterator fam;
  for (fam = families.begin();  fam != families.end();  fam++)
    {
      fam->clearGenotypes();
    }
}


void Pedigree::clearConsistencyVars()
{
  vector<Family>::iterator fam;
  for (fam = families.begin();  fam != families.end();  fam++)
    {
      fam->clearConsistencyVars();
    }
}


// Sets pointers between individuals in the same pedigree
//  (mother, father, and children)
void Pedigree::setFamilyPointers()
{
  // For each family
  vector<Family>::iterator fam;
  for (fam = families.begin();  fam != families.end();  fam++)
    {
      fam->setPedigree(this);
      fam->setIndices();
      fam->setFamilyPointers();
      fam->pedigree_marker_alleles = &marker_alleles;
    }
}


void Family::setPedigree(Pedigree* _pedigree)
{
  pedigree = _pedigree;
}


void Family::setFamilyPointers()
{
  // For each person in the family
  vector<Individual*>::iterator person;
  vector<Individual*>::iterator found;
  for (person = members.begin();  person != members.end();  person++)
    {
      //printf("PERSON: %i, f %i, m %i, s %i\n", person->individual_id, person->father_id, person->mother_id, person->sex);

      // mother
      if ((*person)->mother == NULL && (*person)->mother_id != UNKNOWN)
	{
	  int mother_id = (*person)->mother_id;
	  found = find_individual(members.begin(), 
				  members.end(), mother_id);
	  if (found != members.end())
	    {
	      (*person)->mother = (*found);
	    }
	}

      // father
      if ((*person)->father == NULL && (*person)->father_id != UNKNOWN)
	{
	  int father_id = (*person)->father_id;
	  found = find_individual(members.begin(), 
				  members.end(), father_id);
	  if (found != members.end())
	    {
	      (*person)->father = (*found);
	    }
	}

      // children
      for (found = members.begin();  found != members.end(); found++)
	{
	  //printf("  FOUND: %i, f %i, m %i, s %i\n", (*person)->individual_id, (*person)->father_id, (*person)->mother_id, (*person)->sex);
	  if ((*found)->mother_id == (*person)->individual_id)
	    {
	      if ((*person)->sex == FEMALE)
		{
		  // set mother's child ptr and child's mother pointer
		  (*person)->children.push_back(*found);
		  if ((*found) != NULL){
		    (*found)->mother = (*person);
		  } 
		    
		} else {
		  printf("PED %i ERROR: Sex disagreement between non-female individual %i who is mother of %i\n", family_id, (*person)->individual_id, (*found)->individual_id);
		  exit(-4235);
		}
	    } // end found's mother is person

	  if ((*found)->father_id == (*person)->individual_id)
	    {
	      if ((*person)->sex == MALE)
		{
		  // set father's child ptr and child's father pointer
		  (*person)->children.push_back(*found);
		  if ((*found) != NULL)
		    {
		      (*found)->father = (*person);
		    }
		}
	      else
		{
		  printf("PED %i ERROR: Sex disagreement between non-male individual %i who is father of %i\n", family_id, (*person)->individual_id, (*found)->individual_id);
		  exit(-4236);
		}
	    } // end found's father is person
	}
    }

}





// returns a particular ordering on the members below the founder in the pedigree
// While updating the order vector, the individuals are stored in the same order that they are
// in the family's member list.
void Family::dfs(Individual* parent, int depth)
{

  vector<Individual*>::iterator c;
  for (c = parent->children.begin();  c != parent->children.end();  c++)
    {
      Individual* child = (*c);
      int child_index = child->individual_index;

      // check that order has defined key for this person
      if (order[child_index].individual_index == INVALID)
	{
	  order[child_index].individual_id = child->individual_id;
	  order[child_index].individual_index = child->individual_index;
	  order[child_index].individual = child;
	}
      if (depth > order[child_index].depth)
	{
	  order[child_index].depth = depth;
	  //printf("    child id %i, depth %i\n", child->individual_id, depth);
	  //fflush(stdout);
	}

      dfs(child, depth+1);
    }
}



void Pedigree::setTraversalOrder(bool verbose)
{
  // For each family
  vector<Family>::iterator fam;
  for (fam = families.begin();  fam != families.end();  fam++)
    {
      fam->setTraversalOrder(verbose);
    }
}
 

void Family::setTraversalOrder(bool verbose)
{
  // Allocate the order vector
  for (unsigned int c = 0;  c < members.size();  c++)
    {
      order.push_back(OrderItem(members[c]->individual_id));
    }
  

  // For each founder in the family, determine the top-down order
  // a depth-first search on the pedigree DAG
  for (unsigned int f = 0;  f < founder.size();  f++)
    {
      Individual* founder = this->founder[f];
      
      order[founder->individual_index].individual_index = founder->individual_index;
      order[founder->individual_index].individual_id = founder->individual_id;
      order[founder->individual_index].individual = founder;
      order[founder->individual_index].depth = 0;
      dfs(founder, 1);
    }

  // sort the order vector by decreasing depth
  sort((order).begin(), (order).end(), order_lt());
  
  // set the order in the Individual
  vector<OrderItem>::iterator bottom_up;
  int o = 0;
  for (bottom_up = order.begin();  bottom_up != order.end();  bottom_up++)
  {
    bottom_up->individual->order = o;
    if (verbose){
      printf("      %i] Person %i:%i depth %i\n", o, family_id, bottom_up->individual_id, bottom_up->depth);
      fflush(stdout);
    }
    o++;
  }
}




void Pedigree::setNumMeioses()
{
  vector<Family>::iterator fam;
  for (fam = families.begin();  fam != families.end();  fam++)
    {
      fam->setNumMeioses();
    }
}
void Family::setNumMeioses()
{
  num_meioses = (int) pow(2, 2*nonfounder.size());
}




void Pedigree::setNumAltMeioses()
{
  vector<Family>::iterator fam;
  for (fam = families.begin();  fam != families.end();  fam++)
    {
      fam->setNumAltMeioses();
    }
}

void Family::setNumAltMeioses()
{
  if (nonfounder.empty()){
    num_alt_meioses = 1;
  } else {
    int exponent = 0; // for affected founder
    for (unsigned int i = 0;  i < nonfounder.size();  i++)
      {
	// count meioses that are to unaffecteds
	Individual* indiv = nonfounder[i];
	if (indiv->mother->affection == UNAFFECTED)
	  exponent++;
	if (indiv->father->affection == UNAFFECTED)
	  exponent++;
      }
    //printf("exponent: %i\n", exponent);
    num_alt_meioses = (int) pow(2, exponent);
  }
}






void Pedigree::runFamilyAlgorithm(FamilyAlgorithm& alg)
{
  // For each family
  vector<Family>::iterator fam;
  for (fam = families.begin();  fam != families.end();  fam++)
    {
      alg.run(*fam);
    }
}



void Pedigree::runForwardBlockFamilyAlgorithm(FamilyAlgorithm& alg, int block_size)
{
  int number_of_snps = families.begin()->members[0]->genotypes.size();
  int number_of_blocks = (int) rint(ceil((double) number_of_snps / (double) block_size));

  if (alg.isVerbose()){
    printf("Running forward block family algorithm\n");
    printf("number_of_snps %i, BLOCK_SIZE %i, num of blocks %i\n", number_of_snps, block_size, number_of_blocks);
  }

  alg.setForward();
  for (int block = 0;  block < number_of_blocks;  block++)
    {
      int start_snp = block*block_size;
      int upper_bound_snp = start_snp+block_size;
      if (upper_bound_snp > number_of_snps)
	upper_bound_snp = number_of_snps;
	  

      if (alg.isVerbose())
	printf("\nBLOCK %i, start snp %i, upper_bound snp %i\n", block, start_snp, upper_bound_snp);

      vector<Family>::iterator fam;
      for (fam = families.begin();  fam != families.end();  fam++)
	{
	  alg.setRange(start_snp, upper_bound_snp);
	  alg.run(*fam);
	}

    }
}

void Pedigree::runBackwardBlockFamilyAlgorithm(FamilyAlgorithm& alg, int block_size)
{
  int number_of_snps = families.begin()->members[0]->genotypes.size();
  int number_of_blocks = (int) rint(ceil((double) number_of_snps / (double) block_size));
    


  if (alg.isVerbose()){
    printf("Running backward block family algorithm\n");
    printf("number_of_snps %i, BLOCK_SIZE %i, num of blocks %i\n", number_of_snps, block_size, number_of_blocks);
  }

  alg.setBackward();
  for (int block = number_of_blocks-1;  block >= 0;  block--)
    {
      int start_snp = block*block_size;
      int upper_bound_snp = start_snp+block_size;
      if (upper_bound_snp > number_of_snps)
	upper_bound_snp = number_of_snps;
	  

      if (alg.isVerbose())
	printf("\nBLOCK %i, start snp %i, upper_bound snp %i\n", block, start_snp, upper_bound_snp);

      vector<Family>::iterator fam;
      for (fam = families.begin();  fam != families.end();  fam++)
	{
	  alg.setRange(start_snp, upper_bound_snp);
	  alg.run(*fam);
	}

    }

}



void Pedigree::runForwardBackwardBlockFamilyAlgorithm(FamilyAlgorithm& alg, int block_size)
{
  int number_of_snps = families.begin()->members[0]->genotypes.size();
  int number_of_blocks = (int) rint(ceil((double) number_of_snps / (double) block_size));
    
  if (alg.isVerbose()){
    printf("\n  FORWARD ALG\n-----------------\n");
    printf("number_of_snps %i, BLOCK_SIZE %i, num of blocks %i\n", number_of_snps, block_size, number_of_blocks);
  }



  alg.setForward();
  for (int block = 0;  block < number_of_blocks;  block++)
    {
      int start_snp = block*block_size;
      int upper_bound_snp = start_snp+block_size;
      if (upper_bound_snp > number_of_snps)
	upper_bound_snp = number_of_snps;
	  

      if (alg.isVerbose())
	printf("\nBLOCK %i, start snp %i, upper_bound snp %i\n", block, start_snp, upper_bound_snp);

      vector<Family>::iterator fam;
      for (fam = families.begin();  fam != families.end();  fam++)
	{
	  alg.setRange(start_snp, upper_bound_snp);
	  alg.run(*fam);
	}

    }


  if (alg.isVerbose()){
    printf("\n  BACKWARD ALG\n-----------------\n");
    printf("number_of_snps %i, BLOCK_SIZE %i, num of blocks %i\n", number_of_snps, block_size, number_of_blocks);
  }


  alg.setBackward();
  for (int block = number_of_blocks-1;  block >= 0;  block--)
    {
      int start_snp = block*block_size;
      int upper_bound_snp = start_snp+block_size;
      if (upper_bound_snp > number_of_snps)
	upper_bound_snp = number_of_snps;
	  

      if (alg.isVerbose())
	printf("\nBLOCK %i, start snp %i, upper_bound snp %i\n", block, start_snp, upper_bound_snp);

      vector<Family>::iterator fam;
      for (fam = families.begin();  fam != families.end();  fam++)
	{
	  alg.setRange(start_snp, upper_bound_snp);
	  alg.run(*fam);
	}

    }
}














Individual::Individual()
{
  mother_id = INVALID;
  father_id = INVALID;
  individual_index = INVALID;
  mother = NULL;
  father = NULL;
  any_haplotype = false;
  isSymmetric = true;
  hetero_unknown = false;
  haplotype1 = INVALID;
  haplotype2 = INVALID;
  haplotype1_source = UNKNOWN;
  log_likelihood = 0;
  likelihood=NULL;
  pair_likelihood=NULL;
  sum = 0;
  num_haps = -1;
}


Individual::Individual(int fid, int id, int f, int m, 
		       Gender s, Affection a, long p, int l)
{
  family_id = fid;
  individual_id = id;
  individual_index = INVALID;
  mother_id = m;
  father_id = f;
  mother = NULL;
  father = NULL;
  any_haplotype = false;
  isSymmetric = true;
  hetero_unknown = false;
  sex = s;
  affection = a;
  order = -1;
  file_position = p;
  line_number = l;
  haplotype1 = INVALID;
  haplotype2 = INVALID;
  haplotype1_source = UNKNOWN;
  log_likelihood = 0;
  likelihood=NULL;
  pair_likelihood=NULL;
  sum = 0;
  num_haps = -1;
}

Individual::Individual(const Individual& i)
{
  family_id = i.family_id;
  individual_id = i.individual_id;
  individual_index = i.individual_index;
  father_id = i.father_id;
  mother_id = i.mother_id;
  mother = i.mother;
  father = i.father;
  sex = i.sex;
  affection = i.affection;
  isTyped = i.isTyped;
  order = i.order;
  file_position = i.file_position;
  line_number = i.line_number;
  any_haplotype = i.any_haplotype;
  isSymmetric = i.isSymmetric;
  hetero_unknown = i.hetero_unknown;
  haplotype1 = i.haplotype1;
  haplotype2 = i.haplotype2;
  haplotype1_source = i.haplotype1_source;
  log_likelihood = i.log_likelihood;
  likelihood=NULL;
  pair_likelihood=NULL;
  sum = i.sum;
  num_haps = -1;
  for(unsigned int c = 0;  c < i.genotypes.size(); c++)
    {
      genotypes.push_back(i.genotypes[c]);
    }
  for(unsigned int c = 0;  c < i.inferred_genotypes.size(); c++)
    {
      inferred_genotypes.push_back(i.genotypes[c]);
    }
  for(unsigned int c = 0;  c < i.children.size(); c++)
    {
      children.push_back(i.children[c]);
    }
  for(unsigned int c = 0;  c < i.parental_source.size(); c++)
    {
      parental_source.push_back(i.parental_source[c]);
    }
  for(unsigned int c = 0;  c < i.haplotype_pairs.size(); c++)
    {
      haplotype_pairs.push_back(i.haplotype_pairs[c]);
    }
  for(unsigned int c = 0;  c < i.sampled_haplotype_pairs.size(); c++)
    {
      sampled_haplotype_pairs.push_back(i.sampled_haplotype_pairs[c]);
    }

  num_haps = i.num_haps;
  if (i.likelihood != NULL){
    likelihood = new double[num_haps];
    for(int c = 0;  c < num_haps; c++)
      {
	likelihood[c] = i.likelihood[c];
      }
  }

  if (i.pair_likelihood != NULL){
    pair_likelihood = new double[num_haps*num_haps];
    for(int c = 0;  c < num_haps; c++)
      {
	int row = c*num_haps;
	for (int d = 0;  d < num_haps;  d++)
	  {
	    int index = row + d;
	    pair_likelihood[index] = i.pair_likelihood[index];
	  }
      }
  }
}


const Individual& Individual::operator= (const Individual& i)
{
  family_id = i.family_id;
  individual_id = i.individual_id;
  individual_index = i.individual_index;
  father_id = i.father_id;
  mother_id = i.mother_id;
  mother = i.mother;
  father = i.father;
  sex = i.sex;
  affection = i.affection;
  isTyped = i.isTyped;
  order = i.order;
  file_position = i.file_position;
  line_number = i.line_number;
  any_haplotype = i.any_haplotype;
  isSymmetric = i.isSymmetric;
  hetero_unknown = i.hetero_unknown;
  haplotype1 = i.haplotype1;
  haplotype2 = i.haplotype2;
  haplotype1_source = i.haplotype1_source;
  log_likelihood = i.log_likelihood;
  likelihood=NULL;
  pair_likelihood=NULL;
  sum = i.sum;
  num_haps = -1;
  for(unsigned int c = 0;  c < i.genotypes.size(); c++)
    {
      genotypes.push_back(i.genotypes[c]);
    }
  for(unsigned int c = 0;  c < i.inferred_genotypes.size(); c++)
    {
      inferred_genotypes.push_back(i.genotypes[c]);
    }
  for(unsigned int c = 0;  c < i.children.size(); c++)
    {
      children.push_back(i.children[c]);
    }
  for(unsigned int c = 0;  c < i.parental_source.size(); c++)
    {
      parental_source.push_back(i.parental_source[c]);
    }
  for(unsigned int c = 0;  c < i.haplotype_pairs.size(); c++)
    {
      haplotype_pairs.push_back(i.haplotype_pairs[c]);
    }
  for(unsigned int c = 0;  c < i.sampled_haplotype_pairs.size(); c++)
    {
      sampled_haplotype_pairs.push_back(i.sampled_haplotype_pairs[c]);
    }

  num_haps = i.num_haps;
  if (i.likelihood != NULL){
    likelihood = new double[num_haps];
    for(int c = 0;  c < num_haps; c++)
      {
	likelihood[c] = i.likelihood[c];
      }
  }

  if (i.pair_likelihood != NULL){
    pair_likelihood = new double[num_haps*num_haps];
    for(int c = 0;  c < num_haps; c++)
      {
	int row = c*num_haps;
	for (int d = 0;  d < num_haps;  d++)
	  {
	    int index = row + d;
	    pair_likelihood[index] = i.pair_likelihood[index];
	  }
      }
  }

  return *this;
}







// Returns true if the haplotype is compatible with the 
// known constraints for the given person.
// Returns false otherwise.
bool Individual::haplotypeFits(long haplotype, 
			       const unsigned int start_snp, 
			       const unsigned int upper_bound_snp)
{
  if (haplotype == UNKNOWN_LITERAL) // BK bug fix: Wed 12/19/2007
    return true;                    // BK bug fix: Wed 12/19/2007

  if (haplotype1 == UNKNOWN_LITERAL &&
      haplotype2 == UNKNOWN_LITERAL)
    {
      return true;
    }
  else if (haplotype1 != INVALID && haplotype2 != INVALID)
    {
      // It doesn't work to check for UNKNOWN_LITERAL here,
      // because a child's genotype should match one of the 
      // dictated parent's haplotypes, unless the parent is UU
      //if (haplotype1 == UNKNOWN_LITERAL)
      //  return true;
      //if (haplotype2 == UNKNOWN_LITERAL)
      //  return true;
      

      if (haplotype == haplotype1)
	return true;
      if (haplotype == haplotype2)
	return true;
      return false;
    }


  if (any_haplotype)
    return true;

  bool fit_paternal = true;
  bool fit_maternal = true;




  // check this haplotype against the constraints
  int place = upper_bound_snp - start_snp - 1;
  int value = (long int) pow(2,place);
  for (unsigned int j = start_snp;  j < upper_bound_snp; j++)
    {
      int bin_allele = (int) floor(haplotype / value);
      haplotype = haplotype % value;
      value /= 2;
      int allele = UNKNOWN;
      if (bin_allele == 0)
	allele = FIRST_ALLELE;
      else
	allele = SECOND_ALLELE;



      //if (print_me)
      //printf(" %i (%i,%i ps %i)", allele, genotypes[j].allele1, 
      //       genotypes[j].allele2, parental_source[j]);

      if (parental_source[j] == UNKNOWN || parental_source[j] >= SAME_AS)
	{
	  if (genotypes[j].allele1 != UNKNOWN && genotypes[j].allele2 != UNKNOWN &&
	      genotypes[j].allele1 != allele && genotypes[j].allele2 != allele)
	    {
	      fit_paternal = false;
	      fit_maternal = false;
	      return false;
	    }
	}
      else
	{
	  if (parental_source[j] == MALE)
	    {
	      if (genotypes[j].allele1 != UNKNOWN && 
		  allele != genotypes[j].allele1)
		fit_paternal = false;
	      if (genotypes[j].allele2 != UNKNOWN &&
		  allele != genotypes[j].allele2)
		fit_maternal = false;
	    }
	  else if (parental_source[j] == FEMALE)
	    {
	      if (genotypes[j].allele2 != UNKNOWN &&
		  allele != genotypes[j].allele2)
		fit_paternal = false;
	      if (genotypes[j].allele1 != UNKNOWN &&
		  allele != genotypes[j].allele1)
		fit_maternal = false;
	    }
	}

      if (fit_paternal == false && fit_maternal == false)
	{
	  return false;
	}
    }  

  //if (print_me)
  //printf("\n");


  return true;
}



void Individual::initPairLikelihood(int numofhaps, double init)
{
  num_haps = numofhaps;

  if (pair_likelihood == NULL)
    {
      // Allocate and set
      pair_likelihood = new double[num_haps*num_haps];
      for(unsigned int c = 0;  c < (unsigned) num_haps; c++)
	{
	  for (unsigned int d = 0;  d < (unsigned) num_haps;  d++)
	    {
	      pair_likelihood[c*num_haps+d] = init;
	    }
	}
    }
  else
    {
      // no need to allocate
      for(unsigned int c = 0;  c < (unsigned) num_haps; c++)
	{
	  for (unsigned int d = 0;  d < (unsigned) num_haps;  d++)
	    {
	      pair_likelihood[c*num_haps+d] = init;
	    }
	}
    }
}



////////////////
// returns root locus of SAME_AS set 
int Individual::parental_source_root(int locus)
{
  if (parental_source[locus] < SAME_AS)
    return locus;


  //printf("iterating for locus %i: ", locus);

  // iterate through the pointers until reaching a parent
  int existing_same_index = this->parental_source[locus] - SAME_AS;
  while (this->parental_source[existing_same_index] >= SAME_AS)
    {
      //printf(" %i, ", existing_same_index);
      existing_same_index = this->parental_source[existing_same_index] - SAME_AS;
    }
  //printf(" %i value %i", existing_same_index, parental_source[existing_same_index]);
  assert(existing_same_index >= 0);

  int value = this->parental_source[existing_same_index];
  assert(value == UNKNOWN || value == MALE || value == FEMALE);

  return existing_same_index;
}



////////////////
// returns root locus of SAME_AS set 
void Individual::flatten_parental_source()
{
// #ifdef DEBUG
//   if (family_id == CHECK_F && individual_id == CHECK_I)
//     printf("Flattening... %i:%i\n", family_id, individual_id);
// #endif
  for (unsigned int i = 0;  i < parental_source.size();  i++)
    {
      unsigned int root = parental_source_root(i);
      
#ifdef DEBUG
      //if (family_id == CHECK_F && individual_id == CHECK_I)
      //printf("    SNP %i with value %i has root is %i with value %i\n", i, parental_source[i], root, parental_source[root]);
#endif

      if (parental_source[root] == UNKNOWN && root != i)
	  parental_source[i] = SAME_AS + root;
      else
	  parental_source[i] = parental_source[root];
#ifdef DEBUG
      //if (family_id == CHECK_F && individual_id == CHECK_I)
      //printf("    SNP %i with value %i has root is %i with value %i\n", i, parental_source[i], root, parental_source[root]);
#endif
    }
}






////////////////////////////////
//
//  Checks for consistency between existing parental_source indicators
//  and adds a new source indicator.
//
// precondition: a self-consistent set of parental_source indicators
//               value is MALE or FEMALE
// postcondition: 1) the same set of indicators updated with the newest info
//                2) new set is self-consistent
//                3) program exits if they are any disagreements between
//                   between the new info and the old indicators
//
void Individual::set_value_parental_source(int locus, int new_value)
{
// #ifdef DEBUG
//   if (family_id == CHECK_F && individual_id == CHECK_I)
//     {
//       printf("indiv %i:%i set_value_parental_source(%i,%i)\n", family_id, individual_id, locus, new_value);
//       for (int i = 0;  i < parental_source.size(); i++)
// 	printf(" snp %i, value %i    ", i, parental_source[i]);
//       printf("\n");
//       fflush(stdout);
//     }
// #endif


  if (new_value != MALE  && new_value != FEMALE)
    {
      printf("FATAL ERROR: value (%i) is invalid on call to Individual::set_value_parental_source\n", new_value);
      exit(-5233);
    }

  // check agreement
  int existing_same_index = this->parental_source_root(locus);
  int existing_value = this->parental_source[existing_same_index];
  if (existing_value != UNKNOWN)
    {
      if (existing_value != new_value)
	{
	  printf("FATAL ERROR: Disagreement in %i:%i between parental sources of loci %i with value %i and new value %i\n", family_id, individual_id, existing_same_index, existing_value, new_value);
	  exit(-213);
	}
      this->parental_source[locus] = new_value;
    }
  else // existing_value == UNKNOWN
    {
      this->parental_source[existing_same_index] = new_value;
      this->parental_source[locus] = new_value;
    }

// #ifdef DEBUG
//   if (family_id == CHECK_F && individual_id == CHECK_I)
//     {
//       printf("END indiv %i:%i set_value_parental_source(%i,%i)\n", family_id, individual_id, locus, new_value);
//       for (int i = 0;  i < parental_source.size(); i++)
// 	printf(" snp %i, value %i    ", i, parental_source[i]);
//       printf("\n");
//       fflush(stdout);
//     }
// #endif

}

////////////////////////////////
//
//  Checks for consistency between existing parental_source indicators
//  and adds a new source indicator.
//
// precondition: a self-consistent set of parental_source indicators
//               same_index is [0, #SNPs-1]
// postcondition: 1) the same set of indicators updated with the newest info
//                2) new set is self-consistent
//                3) program exits if they are any disagreements between
//                   between the new info and the old indicators
//
void Individual::set_same_parental_source(int locus, unsigned int same_index)
{
// #ifdef DEBUG
//   if (family_id == CHECK_F && individual_id == CHECK_I)
//   {
//     printf("indiv %i:%i set_same_parental_source(%i,%i)\n", family_id, individual_id, locus, same_index);
//     for (int i = 0;  i < parental_source.size(); i++)
//       printf(" snp %i, value %i    ", i, parental_source[i]);
//     printf("\n");
//     fflush(stdout);
//   }
// #endif


  if (same_index < 0 || same_index > this->genotypes.size())
    {
      printf("FATAL ERROR: same_index (%i) is out of bounds on call to Individual::set_same_parental_source\n", same_index);
      exit(-5233);
    }

  // check agreement
  int existing_same_index = this->parental_source_root(locus);
  int existing_value = this->parental_source[existing_same_index];
  int new_same_index = this->parental_source_root(same_index);
  int new_value = this->parental_source[new_same_index];
  if (existing_value != UNKNOWN && new_value != UNKNOWN)
    {
      //if (family_id == 5 && individual_id == 6) printf("###### -|UNK, -|UNK, newval %i\n", new_value);
      if (existing_value != new_value)
	{
	  printf("\nFATAL ERROR: Disagreement between parental sources of loci %i (%i) and %i (%i)\n", existing_same_index, existing_value, same_index, new_value);
	  exit(-213);
	}
      this->parental_source[locus] = new_value;
    }
  else if (existing_value == UNKNOWN && new_value != UNKNOWN)
    {
      //if (family_id == 1 && individual_id == 7) printf("###### UNK, -|UNK, newval %i\n", new_value);
      this->parental_source[existing_same_index] = new_value;
      this->parental_source[locus] = new_value;
    }
  else if (existing_value != UNKNOWN && new_value == UNKNOWN)
    {
      //if (family_id == 1 && individual_id == 7) printf("###### -|UNK, UNK, newval %i\n", existing_value);
      this->parental_source[new_same_index] = existing_value;
      this->parental_source[locus] = existing_value;
    }
  else
    { // existing_value == UNKNOWN && new_value == UNKNOWN
      int new_parent = existing_same_index;
      if (new_same_index < existing_same_index)
	new_parent = new_same_index;
      //if (family_id == 1 && individual_id == 7) printf("###### UNK, UNK, newparent %i\n", new_parent);

      this->parental_source[new_same_index] = SAME_AS + new_parent;
      this->parental_source[existing_same_index] = SAME_AS + new_parent;
      this->parental_source[locus] = SAME_AS + new_parent;
      this->parental_source[new_parent] = UNKNOWN;
    }

  //if (family_id == 1 && individual_id == 7)
  //{
  //printf("END indiv %i:%i set_same_parental_source(%i,%i)\n", family_id, individual_id, locus, same_index);
  //for (int i = 0;  i < parental_source.size(); i++)
  //  printf(" snp %i, value %i    ", i, parental_source[i]);
  //printf("\n");
  //}
}



















/* void print_block_mendel(Pedigree& pedigree, const char* mendel_file, */
/* 			const unsigned int start_snp, */
/* 			const unsigned int upper_bound_snp) */
/* { */
/*   // print out all inferred haplotypes (along the whole genome) */
/*   FILE* fp = fopen (mendel_file, "w"); */
/*   if (fp == NULL) */
/*     { */
/*       printf("Cannot open (w) %s.", mendel_file); */
/*       exit(-1); */
/*     } */

/*   vector<Family>::iterator fam; */
/*   for (fam = pedigree.families.begin();  fam != pedigree.families.end();  fam++) */
/*     { */
/*       // Print out genotype of every family member */

/*       for (unsigned int i = 0;  i < fam->members.size();  i++) */
/*       { */
/*         Individual* person = fam->members[i]; */
/*         fprintf(fp, "%i %i ", person->family_id, person->individual_id); */
/*         for (unsigned int j = start_snp;  j < upper_bound_snp; j++) */
/*           { */
/*             fprintf(fp, " %i %i [%i]", person->genotypes[j].allele1, person->genotypes[j].allele2, person->parental_source[j]); */
/*           } */
/*         fprintf(fp, "\n"); */
/*       } */
/*     } */

/*   fclose(fp); */
/* } */














//
//  Outputs a linkage file in the Superlink format.
//
void Pedigree::printSuperlinkFile(const char* hap_file)
{
  // print out all inferred haplotypes (along the whole genome)
  FILE* fp = fopen (hap_file, "w");
  if (fp == NULL)
    {
      printf("Cannot open (w) %s.", hap_file);
      exit(-1);
    }


  vector<Family>::iterator fam;
  for (fam = families.begin();  fam != families.end();  fam++)
    {
      vector<Individual*> resorted_members(fam->members.size(), NULL);


      vector<Individual*>::iterator it;
      for (it = fam->members.begin();  it != fam->members.end();  it++)
	{
	  Individual* person = (*it);
	  resorted_members[person->individual_id-1] = person;
	}

      // iterate through the sorted vector
      for (it = resorted_members.begin();  it != resorted_members.end();  it++)
	{
	  Individual* person = (*it);
	  if (person == NULL)
	    {
	      printf("ERROR: null person\n");
	      continue;
	    }

	  int first_child = 0;
	  if (person->children.size() != 0)
	    first_child = person->children[0]->individual_id;
	  int next_sib_father = 0;
	  if (person->father != NULL)
	    {
	      unsigned int index = 0;
	      while (person->father->children[index] != person && index < person->father->children.size())
		{
		  index++;
		}
	      if (index >= person->father->children.size())
		{
		  printf("ERROR: could not find child %i:%i in father's children list\n", person->family_id, person->individual_id);
		  exit(-1);
		}
	      if (index < person->father->children.size()-1)
		next_sib_father = person->father->children[index+1]->individual_id;
	    }
	  int next_sib_mother = 0;
	  if (person->mother != NULL)
	    {
	      unsigned int index = 0;
	      while (person->mother->children[index] != person && index < person->mother->children.size())
		{
		  index++;
		}
	      if (index >= person->mother->children.size())
		{
		  printf("ERROR: could not find child %i:%i in mother's children list\n", person->family_id, person->individual_id);
		  exit(-1);
		}
	      if (index < person->mother->children.size()-1)
		next_sib_mother = person->mother->children[index+1]->individual_id;
	    }


	  fprintf(fp, "%i %i ", person->family_id, person->individual_id);
	  fprintf(fp, "%i %i ", person->father_id, person->mother_id);
	  fprintf(fp, "%i ",    first_child);
	  fprintf(fp, "%i %i ", next_sib_father, next_sib_mother);
	  fprintf(fp, "%i 0 %i ", person->sex, person->affection);
	  fprintf(fp, "%i %i ", person->sex, person->affection);


	  for (unsigned int snp = 0; snp < person->inferred_genotypes.size(); snp++)
	    {
	      long allele1 = person->genotypes[snp].allele1;
	      long allele2 = person->genotypes[snp].allele2;

	      fprintf(fp," %li %li ", allele1, allele2);	      
	    }
	  fprintf(fp, "\n");
	}
    }

  fclose(fp);
}







void Pedigree::printInferredGenotypes(const char* hap_file)
{
  // print out all inferred haplotypes (along the whole genome)
  FILE* fp = fopen (hap_file, "w");
  if (fp == NULL)
    {
      printf("Cannot open (w) %s.", hap_file);
      exit(-1);
    }


  vector<Family>::iterator fam;
  for (fam = families.begin();  fam != families.end();  fam++)
    {
      vector<Individual*>::iterator it;
      for (it = fam->members.begin();  it != fam->members.end();  it++)
	{
	  Individual* person = (*it);
	  fprintf(fp, "%i %i ", person->family_id, person->individual_id);
	  fprintf(fp, "%i %i ", person->father_id, person->mother_id);
	  fprintf(fp, "%i %i ", person->sex, person->affection);
	  
	  for (unsigned int snp = 0; snp < person->inferred_genotypes.size(); snp++)
	    {
	      // If the family is compatible with the no-recomb. hypothesis
	      // we will have values.
	      long allele1 = person->inferred_genotypes[snp].allele1;
	      long allele2 = person->inferred_genotypes[snp].allele2;
	      if (allele1 == UNKNOWN && allele2 == UNKNOWN)
		{
		  // Then the family must have recomb. in order to satisfy 
		  // Mendelian inheritance.  So print out the genotypes.
		  allele1 = person->genotypes[snp].allele1;
		  allele2 = person->genotypes[snp].allele2;
		}
	      
	      fprintf(fp," %li %li ", allele1, allele2);	      
	    }
	  
	  fprintf(fp, "\n");
	}

    }

  fclose(fp);
}








void Pedigree::printInferredGenotypes()
{
  // print out all inferred haplotypes (along the whole genome)
  vector<Family>::iterator fam;
  for (fam = families.begin();  fam != families.end();  fam++)
    {
      vector<Individual*>::iterator it;
      for (it = fam->members.begin();  it != fam->members.end();  it++)
	{
	  Individual* person = (*it);
	  printf("%i %i ", person->family_id, person->individual_id);
	  printf("%i %i ", person->father_id, person->mother_id);
	  printf("%i %i ", person->sex, person->affection);
	  
	  for (unsigned int snp = 0; snp < person->inferred_genotypes.size(); snp++)
	    {
	      // If the family is compatible with the no-recomb. hypothesis
	      // we will have values.
	      long allele1 = person->inferred_genotypes[snp].allele1;
	      long allele2 = person->inferred_genotypes[snp].allele2;
	      if (allele1 == UNKNOWN && allele2 == UNKNOWN)
		{
		  // Then the family must have recomb. in order to satisfy 
		  // Mendelian inheritance.  So print out the genotypes.
		  allele1 = person->genotypes[snp].allele1;
		  allele2 = person->genotypes[snp].allele2;
		}
	      
	      printf(" %li %li ", allele1, allele2);	      
	    }
	  
	  printf("\n");
	}

    }
}



void Pedigree::printGenotypes()
{
  // print out all inferred haplotypes (along the whole genome)
  vector<Family>::iterator fam;
  for (fam = families.begin();  fam != families.end();  fam++)
    {
      vector<Individual*>::iterator it;
      for (it = fam->members.begin();  it != fam->members.end();  it++)
	{
	  Individual* person = (*it);
	  printf("%i %i ", person->family_id, person->individual_id);
	  printf("%i %i ", person->father_id, person->mother_id);
	  printf("%i %i ", person->sex, person->affection);
	  
	  for (unsigned int snp = 0; snp < person->genotypes.size(); snp++)
	    {
	      // If the family is compatible with the no-recomb. hypothesis
	      // we will have values.
	      long allele1 = person->genotypes[snp].allele1;
	      long allele2 = person->genotypes[snp].allele2;
	      
	      printf(" %li %li ", allele1, allele2);	      
	    }
	  
	  printf("\n");
	}

    }
}





void Pedigree::printInferredMultilocusHaps(const char* hap_file, int block_size)
{
  // print out all inferred haplotypes (along the whole genome)
  FILE* fp = fopen (hap_file, "w");
  if (fp == NULL)
    {
      printf("Cannot open (w) %s.", hap_file);
      exit(-1);
    }

  // for each block of 5-snps, print one locus for each person

  int number_of_snps = families.begin()->members[0]->inferred_genotypes.size();
  int number_of_blocks = (int) rint(ceil((double) number_of_snps / (double) block_size));
    
  //printf("number_of_snps %i, BLOCK_SIZE %i, num of blocks %i\n", number_of_snps, BLOCK_SIZE, number_of_blocks);


  vector<Family>::iterator fam;
  for (fam = families.begin();  fam != families.end();  fam++)
    {
      vector<Individual*>::iterator it;
      for (it = fam->members.begin();  it != fam->members.end();  it++)
	{
	  Individual* person = (*it);
	  fprintf(fp, "%i %i ", person->family_id, person->individual_id);
	  fprintf(fp, "%i %i ", person->father_id, person->mother_id);
	  fprintf(fp, "%i %i ", person->sex, person->affection);



	  for (int block = 0;  block < number_of_blocks;  block++)
	    {
	      int start_snp = block*block_size;
	      int upper_bound_snp = start_snp+block_size;
	      if (upper_bound_snp > number_of_snps)
		upper_bound_snp = number_of_snps;

	      assert((unsigned) upper_bound_snp <= person->inferred_genotypes.size());	      


	      //printf("\nBLOCK %i, start snp %i, upper_bound snp %i\n", block, start_snp, upper_bound_snp);

	      // convert each block to an integer haplotype numbered 1...2^size_of_block
	      long int hap1 = 0;
	      long int hap2 = 0;
	      bool all_zero1 = true;
	      bool all_zero2 = true;
	      int last = upper_bound_snp-1;
	      long int bin_place = 1;
	      for (int s = 0; s < upper_bound_snp-start_snp; s++)
		{
		  //printf("%i, ",last-s);
		  long allele1 = person->inferred_genotypes[last-s].allele1;
		  long allele2 = person->inferred_genotypes[last-s].allele2;		 

		  assert (allele1 == SECOND_ALLELE || allele1 == FIRST_ALLELE || allele1 == UNKNOWN);
		  assert (allele2 == SECOND_ALLELE || allele2 == FIRST_ALLELE || allele2 == UNKNOWN);
		  
		  if (allele1 == SECOND_ALLELE)
		    {
		      hap1 += bin_place;
		      all_zero1 = false;
		    }
		  else if (allele1 == FIRST_ALLELE)
		      all_zero1 = false;

		  if (allele2 == SECOND_ALLELE)
		    {
		      hap2 += bin_place;
		      all_zero2 = false;
		    }
		  else if (allele2 == FIRST_ALLELE)
		      all_zero2 = false;

		  
		  bin_place*=2;
		}
	      hap1 += 1;
	      hap2 += 1;
	      if (all_zero1)
		hap1 = 0;
	      if (all_zero2)
		hap2 = 0;

	      //printf("\n");
	      fprintf(fp," %li %li ", hap1, hap2);
	    }
	  fprintf(fp, "\n");
	}
    }

  fclose(fp);
}






void Family::printInconsistencies()
{
  printf("\n");
  printf("Haplotypes for the family %i:\n", family_id);
  for (unsigned int m = 0; m < members.size();  m++)
    {
      Individual* p = members[m];
      printf("   Member %i has haps (%li, %li)", p->individual_id, p->haplotype1, p->haplotype2);
      long int p1 = p->haplotype1;
      long int p2 = p->haplotype2;
      if (p->mother != NULL && p->father != NULL)
	{
	  long int m1 = p->mother->haplotype1;
	  long int m2 = p->mother->haplotype2;
	  long int f1 = p->father->haplotype1;
	  long int f2 = p->father->haplotype2;
	  if ((p1 == m1 || p1 == m2 || m1 == -33 || m2 == -33 || p1 == -33) && (p2 == f1 || p2 == f2 || f1 == -33 || f2 == -33 || p2 == -33))
	    printf("\n");
	  else if ((p2 == m1 || p2 == m2 || m1 == -33 || m2 == -33 || p2 == -33) && (p1 == f1 || p1 == f2|| f1 == -33 || f2 == -33 || p1 == -33))
	    printf("\n");
	  else
	    printf(" ****** inconsistent ******\n");
	}
      else
	printf("\n");
    }
  printf("\n");
}





vector<Individual*>::iterator
 Family::find_individual(vector<Individual*>::iterator curr, 
			 vector<Individual*>::iterator end, 
			 int individual_id)
{
  while (curr != end)
    {
      // check the id for this person
      if ((*curr)->individual_id == individual_id)
	{
	  return (curr);
	}
      curr++;
    }
  return end;
}











// Convert a genotype vector into two integer haplotypes 
// each in the set {0...2^size_of_block-1}.
// The alleles in haplotype_array at positions 
// [start_snp, upper_bound_snp) will be converted into 
// the integers haplotype1 and haplotype2.
void haplotype_to_int(long int& haplotype1, long int& haplotype2, const vector<Marker>& genotype_array, int start_snp, int upper_bound_snp)
{

  long int hap1 = 0;
  long int hap2 = 0;
  bool all_zero1 = true;
  bool all_zero2 = true;
  int last = upper_bound_snp-1;
  long int bin_place = 1;
  for (int s = 0; s < upper_bound_snp-start_snp; s++)
    {
      //printf("%i, ",last-s);
      long allele1 = genotype_array[last-s].allele1;
      long allele2 = genotype_array[last-s].allele2;		 

      assert (allele1 == SECOND_ALLELE || allele1 == FIRST_ALLELE || allele1 == UNKNOWN);
      assert (allele2 == SECOND_ALLELE || allele2 == FIRST_ALLELE || allele2 == UNKNOWN);
		  
      if (allele1 == SECOND_ALLELE)
	{
	  hap1 += bin_place;
	  all_zero1 = false;
	}
      else if (allele1 == FIRST_ALLELE)
	all_zero1 = false;

      if (allele2 == SECOND_ALLELE)
	{
	  hap2 += bin_place;
	  all_zero2 = false;
	}
      else if (allele2 == FIRST_ALLELE)
	all_zero2 = false;

		  
      bin_place*=2;
    }

  if (all_zero1)
    hap1 = UNKNOWN_LITERAL;
  if (all_zero2)
    hap2 = UNKNOWN_LITERAL;

  haplotype1 = hap1;
  haplotype2 = hap2;
}

// Convert a genotype vector into two integer haplotypes 
// each in the set {0...2^size_of_block-1}.
// The alleles in haplotype_array at positions 
// [start_snp, upper_bound_snp) will be converted into 
// the integers haplotype1 and haplotype2.
long int haplotype_to_int(vector<int> allele_array)
{

  int start_snp = 0;
  int upper_bound_snp = allele_array.size();

  long int hap1 = 0;
  bool all_zero1 = true;
  int last = upper_bound_snp-1;
  long int bin_place = 1;
  for (int s = 0; s < upper_bound_snp-start_snp; s++)
    {
      //printf("%i, ",last-s);
      long allele1 = allele_array[last-s];

      assert (allele1 == SECOND_ALLELE || allele1 == FIRST_ALLELE || allele1 == UNKNOWN);
		  
      if (allele1 == SECOND_ALLELE)
	{
	  hap1 += bin_place;
	  all_zero1 = false;
	}
      else if (allele1 == FIRST_ALLELE)
	all_zero1 = false;

      bin_place*=2;
    }

  if (all_zero1)
    hap1 = UNKNOWN_LITERAL;

  return hap1;
}



// The alleles will be written INTO allele_array
void int_to_haplotype(vector<int>& allele_array, const long int haplotype, 
		      const int number_of_snps, 
		      const bool doPrintSTDOUT, FILE* fp)
{
  long int haplotype1 = haplotype;
  int place = number_of_snps-1;
  int value = (long int) pow(2,place);
  for (int j = 0;  j < number_of_snps; j++)
    {
      int bin_allele1 = -1;
      if (haplotype1 != UNKNOWN_LITERAL){
	bin_allele1 = (int) floor(haplotype1 / value);
	haplotype1 = haplotype1 % value;
      }
      
      value /= 2;
      int allele1 = UNKNOWN;
      if (bin_allele1 == 0)
	allele1 = FIRST_ALLELE;
      else if (bin_allele1 == 1)
	allele1 = SECOND_ALLELE;

      allele_array.push_back(allele1);
   
      if (doPrintSTDOUT)
	{
	  printf("%i", allele1);
	}
      if (fp != NULL)
	fprintf(fp," %i ", allele1);
    }
}










double mylog (double x) {
  assert (x >= 0);
  if (x == 0)
    return (-10000);
  else 
    //if (log10(x) < std::numeric_limits<double>::min_exponent10)
    //printf("WARNING: underflow\n");
    return (log(x));
}

int intpow (int x, int y)
{
  int tmp = x;
  int i = 1;

  for (i = 2;  i <= y;  i++)
    tmp = tmp*x;

  return tmp;
}


#endif
