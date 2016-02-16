//////////////////////////////////////////////////
//
// Converts from the typical linkage format into that required by Superlink.
//
// Superlink format requires more pointers describing the pedigree graph:
//     1) a pedigree number
//     2) an individual identification number, or id
//     3) father's id number
//     4) mother's id number
//     5) first child
//     6) next sibling by father
//     7) next sibling by mother)
//     8) sex (male 1, female 2)  -- 72 male, 60 female
//     9) zero
//     10) affection status (affected 2, unaffected 1, unknown 0)
//
//////////////////////////////////////////////////




#include <vector>
#include <algorithm>
#include <string>
#include <set>
#include <queue>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <time.h>



using namespace std;




#define FALSE 0
#define TRUE 1


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


#define DOUBLEDIGIT -11
#define INVALID_ALLELE -22
#define OUTOFBOUNDS -33

#define CHECK_F 5
#define CHECK_I 6

//#define DEBUG 1  // comment out to remove debugging

#define USE_UNKNOWN_LITERAL 1


int verbose = false;
int block_size = -1;
const double log_half = log(0.5);


typedef int Gender;
typedef int Affection;




class Marker;
class Family;
class Individual;

vector<Family> pedigrees;








class Marker
{
 public:
  Marker(){
    allele1 = UNKNOWN;
    allele2 = UNKNOWN;
  }

  // make sure alleles are sorted when storing them as markers
  Marker(int a, int b){
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
  int insertAllele(int a)
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

  int allele1;
  int allele2;
};






class Individual
{
public:
  Individual(){
    mother_id = INVALID;
    father_id = INVALID;
    individual_index = INVALID;
    mother = NULL;
    father = NULL;
    any_haplotype = FALSE;
    isSymmetric = TRUE;
    haplotype1 = INVALID;
    haplotype2 = INVALID;
    haplotype1_source = UNKNOWN;
    log_likelihood = 0;
  }
  Individual(int fid, int id, int f, int m, Gender s, Affection a, long p, int l){
    family_id = fid;
    individual_id = id;
    individual_index = INVALID;
    mother_id = m;
    father_id = f;
    mother = NULL;
    father = NULL;
    any_haplotype = FALSE;
    isSymmetric = TRUE;
    sex = s;
    affection = a;
    file_position = p;
    line_number = l;
    haplotype1 = INVALID;
    haplotype2 = INVALID;
    haplotype1_source = UNKNOWN;
    log_likelihood = 0;
  }
  Individual(const Individual& i){
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
    file_position = i.file_position;
    line_number = i.line_number;
    disease_marker = disease_marker;
    any_haplotype = i.any_haplotype;
    isSymmetric = i.isSymmetric;
    haplotype1 = i.haplotype1;
    haplotype2 = i.haplotype2;
    haplotype1_source = i.haplotype1_source;
    log_likelihood = i.log_likelihood;
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
    for(unsigned int c = 0;  c < i.likelihood.size(); c++)
      {
	likelihood.push_back(i.likelihood[c]);
      }
    for(unsigned int c = 0;  c < i.pair_likelihood.size(); c++)
      {
	pair_likelihood.push_back(i.pair_likelihood[c]);
      }
  }
  const Individual& operator= (const Individual& i){
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
    file_position = i.file_position;
    line_number = i.line_number;
    disease_marker = disease_marker;
    any_haplotype = i.any_haplotype;
    isSymmetric = i.isSymmetric;
    haplotype1 = i.haplotype1;
    haplotype2 = i.haplotype2;
    haplotype1_source = i.haplotype1_source;
    log_likelihood = i.log_likelihood;
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
    for(unsigned int c = 0;  c < i.likelihood.size(); c++)
      {
	likelihood.push_back(i.likelihood[c]);
      }
    for(unsigned int c = 0;  c < i.pair_likelihood.size(); c++)
      {
	pair_likelihood.push_back(i.pair_likelihood[c]);
      }
    return *this;
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

  Marker disease_marker;


  // one source for each marker
  // A glorified set data structure
  // 0 -> unknown
  // 1 -> allele a (a < b) from father
  // 2 -> allele a (a < b) from mother
  // SAME_AS -> parental_source[j] = parental_source[offset]
  vector<int> parental_source; 
  int isSymmetric;

  // Entries come in pairs.
  // All haplotypes must be valid [0..2^S-1] 
  // where S is the number of SNPs.
  // The first haplotype from the pair is the 
  // one from the father.
  vector<long> haplotype_pairs;


  long int haplotype1;
  long int haplotype2;
  int haplotype1_source; // MALE, FEMALE, or UNKNOWN

  double log_likelihood;


  // allocated only in people who are  not founders or children of founders
  vector<double> likelihood; // normalized probabilities
  // allocated only in children of founders
  vector<vector<double> > pair_likelihood;  // normalized probabilities

};







class Family
{
public:
  Family(int id){ 
    family_id = id; 
    isCompatible = FALSE;
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
    


  vector<Individual*>::iterator find_individual(vector<Individual*>::iterator curr, 
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



  int family_id;
  vector<Individual*> members;

  vector<Individual*> founder; // index of the child in the vector of individuals
  vector<Individual*> nonfounder; // index of the child in the vector of individuals
  vector<Individual*> typed; // index of the child in the vector of individuals
  vector<Individual*> untyped; // index of the child in the vector of individuals


  set<long> candidate_haplotypes;

  //vector<OrderItem> order;  // bottom up order
  //vector<FounderPairDLists*> descent_lists;

  bool isCompatible;
};











int read_allele (char* line, int length, int& index)
{
  while (index < length && isspace(*(line + index)) ){
    index++;
  }
  //printf("read_allele(len %i, idx %i) line[idx] = %c\n", length, index, *(line + index));

  if (index >= length){ 
    //printf("Returning OUTOFBOUNDS due to %i >= %i", index, length);
    return OUTOFBOUNDS;
  } 

  string a;
  a = a + *(line+index);
  int len = 1;
  while (!isspace(*(line+index+len)) && *(line+index+len) != '\0'){
    a = a + *(line+index+len);
    len++;
  }
  index+=len;



  int allele = atoi(a.c_str());


  // For valid SNP alleles
  if (len > 1){
    printf("double digit allele: |%s|\n", a.c_str());
    return DOUBLEDIGIT;
  }
  // check for valid SNP allele
  //if (allele != UNKNOWN && allele != FIRST_ALLELE && allele != SECOND_ALLELE){
  //  return INVALID_ALLELE;
  //}


  // check for arbitrary marker allele (i.e. non-SNP markers)
  if (allele < UNKNOWN){
    return INVALID_ALLELE;
  }

  while (index < length && isspace(*(line + index)) ){
    index++;
  }

  return allele;
}





void read_pedigree_file(const char* in_file)
{


  // read the graph  
   FILE* fp = fopen (in_file, "r");
   if (fp == NULL)
     {
      printf("\nERROR: Cannot open %s.\n\n", in_file);
      exit(-1);
   }



  char line [1024];



  char column[256];
  int line_number = 1;
  long file_position = ftell(fp);
  while ( fscanf (fp, "%s ", column) != EOF)
      {
	//printf("%i]", line_number);
	int famID = atoi(column);
	//printf(" %i", famID);

	fscanf (fp, "%s ", column);
	int indivID = atoi(column);
	//printf(" %i", indivID);

	fscanf (fp, "%s ", column);
	int faID = atoi(column);
	//printf(" %i", faID);

	fscanf (fp, "%s ", column);
	int moID = atoi(column);
	//printf(" %i", moID);

	int isFounder = FALSE;
	if (faID == UNKNOWN && moID == UNKNOWN)
	  {
	    isFounder = TRUE;
	  }


	fscanf (fp, "%s ", column);
	int sex = atoi(column);
	//printf(" %i", sex);

	fscanf (fp, "%s", column);
	int aff = atoi(column);
	//printf(" %i", aff);
	//printf("\n");
	//fflush(stdout);



	// Find the family or make it
	Family newfam(famID);
	vector<Family>::iterator fam;
	fam = find(pedigrees.begin(), pedigrees.end(), newfam);
	if (fam == pedigrees.end()){
	  // make family
	  pedigrees.push_back(newfam);
	  fam = pedigrees.end() - 1;
	}

	//printf(" newfam  ");
	//fflush(stdout);


	// insert individual into pedigree
	Individual* newperson = new Individual(famID, indivID, faID, moID, sex, aff, file_position, line_number);
	vector<Individual*>::iterator person;
	(fam->members).push_back(newperson);
	person = fam->members.end() - 1;
	if (isFounder){
	  (fam->founder).push_back(newperson);
	} else {
	  (fam->nonfounder).push_back(newperson);
	}




	//printf(" newper  ");
	//fflush(stdout);


	// Read all the alleles
	int isTyped = FALSE;
	
	
	char* notNull = fgets(line, 1024, fp);
	int length = strlen(line);
	int odd_allele = INVALID_ALLELE;


	// in case it didn't read an end-of-line character
	while (length > 0 && notNull != NULL)
	  {
	    // for each SNP
	    int index = 0;
	    while (index  < length) // index++ occurs in read_allele
	      {
		int allele1 = INVALID_ALLELE;
		if (odd_allele != INVALID_ALLELE){
		  allele1 = odd_allele;
		  odd_allele = INVALID_ALLELE;
		} else {
		  //printf("read_allele(len %i, idx %i) line[idx] = |%c|\n", length, index, line[index]);
		  allele1 = read_allele(&line[0], length, index);
		}
		if (index >= length && allele1 != OUTOFBOUNDS){
		  //printf("Odd allele: allele1 %i\n", allele1);
		  odd_allele = allele1;
		  break;
		}
		int allele2 = read_allele(&line[0], length, index);


		if (allele1 == DOUBLEDIGIT || allele2 == DOUBLEDIGIT){
		  printf(" INPUT ERROR: Double digit allele read on line %i\n", line_number);
		  exit(-1);
		}
		if (allele1 == INVALID_ALLELE || allele2 == INVALID_ALLELE){
		  printf(" INPUT ERROR: invalid allele (not %i or %i) read on line %i\n", FIRST_ALLELE, SECOND_ALLELE, line_number);
		  exit(-1);
		}
		if ((allele1 == OUTOFBOUNDS && allele2 != OUTOFBOUNDS)
		     || (allele2 == OUTOFBOUNDS && allele1 != OUTOFBOUNDS)){
		  printf(" FATAL ERROR: odd number of alleles in fgets() on line %i\n", line_number);
		  printf("   allele1 %i   allele2 %i\n", allele1, allele2);
		  printf("   at locus # %i\n", ((*person)->genotypes.size()+1) );
		  printf("   index %i, length %i\n", index, length);
		  printf("   odd_allele %i\n", odd_allele);
		  exit(-1);
		}

		if (allele1 != OUTOFBOUNDS && allele2 != OUTOFBOUNDS)
		  {
		    if (!isTyped && (allele1 != UNKNOWN || allele2 != UNKNOWN))
		      {
			isTyped = TRUE;
		      }

		    //printf("  %i/%i", allele1, allele2); fflush(stdout);
		    Marker snp(allele1, allele2);
		    (*person)->genotypes.push_back(snp);
		    (*person)->inferred_genotypes.push_back(Marker(UNKNOWN, UNKNOWN));
		    (*person)->parental_source.push_back(UNKNOWN);
		  }
	      }


	    // read the next part of the buffer
	    if (line[length-1] != '\n' && notNull != NULL){
	      //printf("fgets\n");
	      notNull = fgets(line, 1024, fp);
	      length = strlen(line);
	    } else {
	      length = 0;
	    }
	  }

	//printf("read %i:%i\n",  (*person)->family_id, (*person)->individual_id);
	//printf("  typed? %i\n", isTyped); fflush(stdout);
	//printf("  num snps? %i\n", (*person)->genotypes.size()); fflush(stdout);




	// add individual to typed and untyped arrays
	if (isTyped){
	  (fam->typed).push_back( *person );
	} else {
	  (fam->untyped).push_back( *person );
	}
	(*person)->isTyped = isTyped;





	line_number++;

	if (feof(fp) != 0)
	  break;


	file_position = ftell(fp);	
      } // end while not EOF



  fclose(fp); // close in file 

}











// Sets pointers between individuals in the same pedigree
//  (mother, father, and children)
void set_pedigree_pointers()
{
  // For each family
  vector<Family>::iterator fam;
  for (fam = pedigrees.begin();  fam != pedigrees.end();  fam++)
    {
      fam->setIndices();

      // For each person in the family
      vector<Individual*>::iterator person;
      vector<Individual*>::iterator found;
      for (person = fam->members.begin();  person != fam->members.end();  person++)
	{
	  //printf("PERSON: %i, f %i, m %i, s %i\n", person->individual_id, person->father_id, person->mother_id, person->sex);

	  // mother
	  if ((*person)->mother == NULL && (*person)->mother_id != UNKNOWN)
	    {
	      int mother_id = (*person)->mother_id;
	      found = fam->find_individual(fam->members.begin(), 
					   fam->members.end(), mother_id);
	      if (found != fam->members.end())
		{
		  (*person)->mother = (*found);
		}
	    }

	  // father
	  if ((*person)->father == NULL && (*person)->father_id != UNKNOWN)
	    {
	      int father_id = (*person)->father_id;
	      found = fam->find_individual(fam->members.begin(), 
					   fam->members.end(), father_id);
	      if (found != fam->members.end())
		{
		  (*person)->father = (*found);
		}
	    }

	  // children
	  for (found = fam->members.begin();  found != fam->members.end(); found++)
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
		      printf("PED %i ERROR: Sex disagreement between non-female individual %i who is mother of %i\n", fam->family_id, (*person)->individual_id, (*found)->individual_id);
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
		      printf("PED %i ERROR: Sex disagreement between non-male individual %i who is father of %i\n", fam->family_id, (*person)->individual_id, (*found)->individual_id);
		    }
		} // end found's father is person
	    }
	}
    }


}







// Superlink format requires more pointers describing the pedigree graph:
//     1) a pedigree number
//     2) an individual identification number, or id
//     3) father's id number
//     4) mother's id number
//     5) first child
//     6) next sibling by father
//     7) next sibling by mother
//     8) sex (male 1, female 2)  -- 72 male, 60 female
//     9) zero
//     10) affection status (affected 2, unaffected 1, unknown 0)

void print_superlink_file(const char* hap_file)
{
  // print out all inferred haplotypes (along the whole genome)
  FILE* fp = fopen (hap_file, "w");
  if (fp == NULL)
    {
      printf("Cannot open (w) %s.", hap_file);
      exit(-1);
    }


  vector<Family>::iterator fam;
  for (fam = pedigrees.begin();  fam != pedigrees.end();  fam++)
    {
      vector<Individual*>::iterator it;
      for (it = fam->members.begin();  it != fam->members.end();  it++)
	{
	  Individual* person = (*it);


	  int first_child = 0;
	  if (person->children.size() != 0)
	    first_child = person->children[0]->individual_id;
	  int next_sib_father = 0;
	  if (person->father != NULL)
	    {
	      int index = 0;
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
	      int index = 0;
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


	  for (unsigned int snp = 0; snp < person->inferred_genotypes.size(); snp++)
	    {
	      long allele1 = person->genotypes[snp].allele1;
	      long allele2 = person->genotypes[snp].allele2;

	      fprintf(fp," %i %i ", allele1, allele2);	      
	    }
	  fprintf(fp, "\n");
	}
    }

  fclose(fp);
}






void deleteped()
{
  // For each family
  vector<Family>::iterator fam;
  for (fam = pedigrees.begin();  fam != pedigrees.end();  fam++)
    {
      // For each person in the family
      vector<Individual*>::iterator person;
      for (person = fam->members.begin();  person != fam->members.end();  person++)
	{
	  delete (*person);
	}


      // Delete the descent_lists
      //vector<FounderPairDLists*>::iterator fdlist;
      //for(fdlist = fam->descent_lists.begin();  fdlist != fam->descent_lists.end();  fdlist++)
      //{
      //  delete (*fdlist);
      //}
    }  
}




//////////////////////////////
//
//
int main (int argc, char *argv[])
{
  // read arguements
  if (argc < 3)  // argc is # argvs [1..n] with prog as argv[0]
    {
      printf ("\nUSAGE: %s [in pedigree file] [out pedigree file]\n\n", argv[0]);
      printf ("[in pedigree file]   INPUT pedigree (genotypes ignored) in linkage format:\n");
      printf ("            FamID IndID FaID MoID Sex Aff [LOCI ...]\n");
      printf ("[out pedigree file]  OUT pedigree file name with SUPERLINK format.\n");
      printf ("\n");

      exit(-1);
    }

  
  const char* ped_file = argv[1];
  const char* out_ped_file = argv[2];



  printf("\nPARAMS:  in ped = %s,\n", ped_file);
  printf("         out ped = %s\n\n", out_ped_file);





  printf("READING input pedigree...\n");
  read_pedigree_file(ped_file);
  set_pedigree_pointers();
  //printf("Finished reading file\n");
  //fflush(stdout);


  printf("PRINTING SUPERLINK FORMAT...\n");


  print_superlink_file(out_ped_file);




  deleteped();
}
