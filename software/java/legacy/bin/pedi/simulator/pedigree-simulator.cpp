
///////////////////////////////////////////
//
//  Given a pedigree of individuals, population haplotypes, and 
//  recombination rates, the simulator draws haplotypes for the 
//  founders and simulate inheritance for their descendents.
//
//  If a disease model is given, a disease SNP is chosen, 
//  and the affection status is simulated based on that SNP.
//  Otherwise, the disease status is unaltered.
//
///////////////////////////////////////////


#include <vector>
#include <algorithm>
#include <string>
#include <set>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <time.h>
#include <string.h>


using namespace std;


#define CASE_FREQUENCY 0.2  // for rejection sampling


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


//#define DEBUG 1  // comment out to remove debugging







typedef int Gender;
typedef int Affection;



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
  }
  Individual(int fid, int id, int f, int m, Gender s, Affection a, long p, int l){
    family_id = fid;
    individual_id = id;
    individual_index = INVALID;
    mother_id = m;
    father_id = f;
    mother = NULL;
    father = NULL;
    sex = s;
    affection = a;
    file_position = p;
    line_number = l;
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
    for(unsigned int c = 0;  c < i.genotypes.size(); c++)
      {
	genotypes.push_back(i.genotypes[c]);
      }
    for(unsigned int c = 0;  c < i.children.size(); c++)
      {
	children.push_back(i.children[c]);
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
    for(unsigned int c = 0;  c < i.genotypes.size(); c++)
      {
	genotypes.push_back(i.genotypes[c]);
      }
    for(unsigned int c = 0;  c < i.children.size(); c++)
      {
	children.push_back(i.children[c]);
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


  vector<Individual*> children; 
  vector<Marker> genotypes;
  vector<Marker> identity;  // name of ancestor alleles that these were inherited from (i.e. IBD state)

  Marker disease_marker;
  Marker disease_identity; // disease SNP IBD state

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


class Family
{
public:
  Family(int id){ family_id = id; }

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


  vector<OrderItem> order;  // bottom up order

};



vector<Family> pedigrees;
vector<vector<int> > founder_haplotypes;
vector<double> recombination_rates;





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
  if (allele != UNKNOWN && allele != FIRST_ALLELE && allele != SECOND_ALLELE){
    return INVALID_ALLELE;
  }


  // check for arbitrary marker allele (i.e. non-SNP markers)
  //if (allele < UNKNOWN){
  //  return INVALID_ALLELE;
  //}

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


		    // We don't want to read the data, so do nothing
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
	//printf("  typed? %i\n", isTyped); fflush(stdout);



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



void read_recombination_file(const char* in_file, int number_of_snps)
{


  // read the graph  
   FILE* fp = fopen (in_file, "r");
   if (fp == NULL)
     {
      printf("\nERROR: Cannot open %s.\n\n", in_file);
      exit(-1);
   }


   char column[256];
   double rec_rate = -1.0;
   int count_rates = 0;
   //printf("Recombination rates:  ");
   while ( fscanf (fp, "%s ", column) != EOF)
     {
       rec_rate = atof(column);
       recombination_rates.push_back(rec_rate);
       count_rates++;
       //printf(" %f", rec_rate);
     }
   //printf("\n");

   if (count_rates > number_of_snps-1)
     {
       printf("ERROR: number of recomination rates exceeds one less than the number of loci\n");
       exit(-534);
     }
   if (count_rates < number_of_snps-1)
     {
       printf("ERROR: number of recomination rates is fewer than one less than the number of loci\n");
       exit(-534);
     }
   
  fclose(fp);
}



int read_founder_file(const char* in_file)
{
  // read the graph  
   FILE* fp = fopen (in_file, "r");
   if (fp == NULL)
     {
      printf("\n ERROR: Cannot open %s.\n\n", in_file);
      exit(-1);
   }



  char line [1024];
  int line_number = 1;
  int number_of_loci = 0;


  while (feof(fp) == 0)
    {
      
      char* notNull = fgets(line, 1024, fp);
      int length = strlen(line);



      int locus_count = 0;
      vector<int> haplotype;
      
      // in case it didn't read an end-of-line character
      while (length > 0 && notNull != NULL)
	{
	  // for each SNP
	  int index = 0;
	  while (index  < length) // index++ occurs in read_allele
	    {
	      int allele1 = read_allele(&line[0], length, index);

	      //printf("index %i, length %i, notNull %i, allele %i \n", index, length, (notNull!=NULL), allele1);


	      if (allele1 == DOUBLEDIGIT){
		printf(" INPUT ERROR: Double digit allele read on line %i\n", line_number);
		exit(-1);
	      }
	      if (allele1 == INVALID_ALLELE){
		printf(" INPUT ERROR: invalid allele (not %i or %i) read on line %i\n", FIRST_ALLELE, SECOND_ALLELE, line_number);
		exit(-1);
	      }
	      if (allele1 != OUTOFBOUNDS && allele1 != 0 && allele1 != 1)
		{
		  printf(" INPUT ERROR: invalid allele %i (not %i or %i) read on line %i\n", allele1, 0, 1, line_number);
		  exit(-1);
		}

	      if (allele1 != OUTOFBOUNDS)
		{
		  locus_count++;
		  haplotype.push_back(allele1);
		  //printf("  %i", allele1); fflush(stdout);
		}

	    } // end while not end of buffer


	  // read the next part of the buffer
	  if (line[length-1] != '\n' && notNull != NULL){
	    //printf("fgets\n");
	    notNull = fgets(line, 1024, fp);
	    length = strlen(line);
	  } else {
	    length = 0;
	  }

	  //printf("index %i, length %i, notNull %i \n", index, length, (notNull!=NULL));

	} // for each line


      if (line_number == 1)
	  number_of_loci = locus_count;

      
      if (locus_count != number_of_loci && locus_count != 0)
	{
	  printf("\nERROR: fewer loci (%i) in line %i than in line 1 (expected %i)\n", locus_count, line_number, number_of_loci);
	  exit(-1);
	}
      else if (locus_count != 0)
	{
	  // push this haplotype into a list of haplotypes
	  founder_haplotypes.push_back(haplotype);
	  //printf("   %i]", (founder_haplotypes.size()+1));
	  //for (int a = 0;  a < haplotype.size();  a++)
	  //    printf(" %i", haplotype[a]);
	  //printf("\n");
	}

      
      line_number++;

    } // while not end of file

      
  fclose(fp); // close in file 


  return number_of_loci;
};



//////////////////
// snp must be a number in range [0, m-1]
//   where m is the number of loci
//
// post-condition: sets minor_allele_freq and minor allele
//                 for snp
void get_founder_allele_freq(unsigned int snp, double& minor_allele_freq, int& minor_allele)
{
  if (snp < 0 || snp > founder_haplotypes[0].size()-1)
    {
      printf("ERROR: snp %i out of range [0,%i]\n", snp, (founder_haplotypes[0].size()-1));
      exit(-1);
    }

  double count = 0.0;
  double total = 0.0;
  for (unsigned int i = 0;  i < founder_haplotypes.size();  i++)
    {
      count += founder_haplotypes[i][snp];
      total += 1;
    }

  double f = count / total;
  printf("Disease allele freq: %f / %f = %f\n", count, total, f);

  double allele_one_freq = count / total;
  if (allele_one_freq > 0.5)
    {
      minor_allele_freq = 1-allele_one_freq;
      minor_allele = 0;
    }
  else
    {
      minor_allele_freq = allele_one_freq;
      minor_allele = 1;
    }
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




void allocate_markers(int number_of_snps)
{
  // For each family
  vector<Family>::iterator fam;
  for (fam = pedigrees.begin();  fam != pedigrees.end();  fam++)
    {
      //printf("FAMILY %i\n", fam->family_id);

      // For each person in the family
      vector<Individual*>::iterator person;
      for (person = fam->members.begin();  person != fam->members.end();  person++)
	{
	  // Allocate memory for all the markers we will simulate
	  for (int marker = 0;  marker < number_of_snps;  marker++)
	    {
	      Marker snp(UNKNOWN, UNKNOWN);
	      (*person)->genotypes.push_back(snp);
	      (*person)->identity.push_back(snp);
	    }
	}
    }
}












// returns a particular ordering on the members below the founder in the pedigree
// While updating the order vector, the individuals are stored in the same order that they are
// in the family's member list.
int pedigree_dfs(vector<Family>::iterator fam, 
		  Individual* parent, 
		  int depth)
{
  int min_depth = 999999999;
  vector<Individual*>::iterator c;
  for (c = parent->children.begin();  c != parent->children.end();  c++)
    {
      Individual* child = (*c);
      int child_index = child->individual_index;

      // check that fam->order has defined key for this person
      if (fam->order[child_index].individual_id == UNKNOWN)
	{
	  fam->order[child_index].individual_id = child->individual_id;
	  fam->order[child_index].individual_index = child->individual_index;
	  fam->order[child_index].individual = child;

	  fam->order[child_index].depth = pedigree_dfs(fam, child, depth+1);
	}
      if (fam->order[child_index].depth < min_depth)
	{
	  min_depth = fam->order[child_index].depth;
	}
    }

  return min_depth-1;

}


void set_traversal_orders()
{
  // For each family
  vector<Family>::iterator fam;
  for (fam = pedigrees.begin();  fam != pedigrees.end();  fam++)
    {
      //printf("FAMILY %i\n", fam->family_id);

      // Allocate the order vector
      for (unsigned int c = 0;  c < fam->members.size();  c++)
	{
	  fam->order.push_back(OrderItem());
	}

      int max_depth = fam->members.size();
      for (unsigned int c = 0; c < fam->members.size(); c++)
	{
	  Individual* leaf = fam->members[c];
	  if (leaf->children.size() == 0)
	    {
	      fam->order[leaf->individual_index].individual_index = leaf->individual_index;
	      fam->order[leaf->individual_index].individual_id = leaf->individual_id;
	      fam->order[leaf->individual_index].individual = leaf;
	      fam->order[leaf->individual_index].depth = max_depth;
	    }
	}


      // For each founder in the family, determine the top-down order
      // a depth-first search on the pedigree DAG
      for (unsigned int f = 0;  f < fam->founder.size();  f++)
	{
	  Individual* founder = fam->founder[f];
	  fam->order[founder->individual_index].individual_index = founder->individual_index;
	  fam->order[founder->individual_index].individual_id = founder->individual_id;
	  fam->order[founder->individual_index].individual = founder;

	  //printf("ped dfs start founder %i:%i....\n", founder->family_id, founder->individual_id);
	  fam->order[founder->individual_index].depth = pedigree_dfs(fam, founder, 1);
	  //printf("ped dfs finished.\n");
	}

      // sort the order vector by decreasing depth
      sort((fam->order).begin(), (fam->order).end(), order_lt());


      //vector<OrderItem>::iterator bottom_up;
      //for (bottom_up = fam->order.begin();  bottom_up != fam->order.end();  bottom_up++)
      //{
      //  printf("      Person %i:%i depth %i\n", fam->family_id, bottom_up->individual_id, bottom_up->depth);
      //  fflush(stdout);
      //}
    }
}




















int sample_phenotypes(int disease_locus, double p, int minor_allele, 
		      double f0, double f1, double f2,
		      double fraction_cases)
{
  printf("sample_pheno(dl %i, p %f, ma %i, f0 %f, f1 %f, f2 %f, fc %f)\n", disease_locus, p, minor_allele, f0, f1, f2, fraction_cases);


  int number_of_cases = 0;
  int number_of_individuals = 0;

  // For each family
  vector<Family>::iterator fam;
  for (fam = pedigrees.begin();  fam != pedigrees.end();  fam++)
    {
      int num_indiv = fam->members.size();
      number_of_individuals += num_indiv;
      int required_cases = (int) floor(fraction_cases * num_indiv);

      // setting a non-zero value is importance sampling, so the probability of the sample must be asjusted
      double required_founder_alleles = 0;  //fraction_cases*fam->founder.size()*2*1.2;


      printf("FAMILY %i, with %i required cases, founder disease alleles %f\n", fam->family_id, required_cases, required_founder_alleles);



      // REJECTION SAMPLING
      int isRejected = true;
      while (isRejected)
	{
	  int fam_num_cases = 0;
	  int founder_num_cases = 0;
	  int founder_num_disease_allele = 0;

#ifdef DEBUG
	  printf("    ... Founders\n");
#endif

	  //while (founder_num_disease_allele == 0)
	  while (founder_num_disease_allele <= required_founder_alleles)
	    {
	      vector<Individual*>::iterator person;
	      founder_num_disease_allele = 0;
	      founder_num_cases = 0;
	      int identity_counter=1;
	      for (person = fam->founder.begin();  person != fam->founder.end();  person++)
		{
		  //printf("  (sample phenotypes) Processing indiv %i:%i \n", fam->family_id, (*person)->individual_id);
		  //fflush(stdout);
		  
		  Individual* founder = (*person);
		  
		  int allele1 = UNKNOWN;
		  int allele2 = UNKNOWN;
		  double r = 0.0;
		  

		      
		  allele1 = 1-minor_allele;
		  allele2 = 1-minor_allele;
		  r = drand48();
		  if (r < p)
		    {
		      allele1 = minor_allele;
		      founder_num_disease_allele++;
		    }
		  r = drand48();
		  if (r < p)
		    {
		      allele2 = minor_allele;
		      founder_num_disease_allele++;
		    }
		  
		  double prob_disease;
		  if (allele1 == minor_allele && allele2 == minor_allele)
		    prob_disease = f2;
		  else if (allele1 != minor_allele && allele2 != minor_allele)
		    prob_disease = f0;
		  else
		    prob_disease = f1;

		  r = drand48();		  
		  founder->affection = UNAFFECTED;
		  if (r < prob_disease)
		    {
		      founder->affection = AFFECTED;
		      founder_num_cases++;
		    }
		  
		  founder->disease_marker.allele1 = allele1+1;
		  founder->disease_marker.allele2 = allele2+1;
		  founder->disease_identity.allele1 = identity_counter;
		  founder->disease_identity.allele2 = identity_counter+1;
		  identity_counter+=2;
		  //printf("         founder %i] disease alleles (%i,%i) phenotype %i\n", founder->individual_id, founder->disease_marker.allele1, founder->disease_marker.allele2, founder->affection);	  
		}

	    } // end while no founder has the disease allele

	  fam_num_cases += founder_num_cases;

#ifdef DEBUG
	  printf("    ... NON-founders\n");
#endif

	  // for non-founder individual (in top-down order)
	  vector<OrderItem>::iterator top_down;
	  for (top_down = fam->order.end()-1;  top_down != fam->order.begin()-1;  top_down--)
	    {
	      //printf("  (sample phenotypes) Processing indiv %i:%i \n", fam->family_id, top_down->individual_id);
	      //fflush(stdout);
		  
		  
	      Individual* individual = top_down->individual;

	      int allele1 = UNKNOWN;
	      int allele2 = UNKNOWN;
	      int id1 = UNKNOWN;
	      int id2 = UNKNOWN;
	      double r = 0.0;

	      // is not founder
	      if (individual->father_id != UNKNOWN && individual->mother_id != UNKNOWN)
		{

		  if (individual->father == NULL || individual->mother == NULL)
		    {
		      printf("ERROR: a parent pointer is NULL\n");
		      exit(-24);
		    }
		  
		  allele1 = individual->father->disease_marker.allele1;
		  id1 = individual->father->disease_identity.allele1;
		  allele2 = individual->mother->disease_marker.allele1;
		  id2 = individual->mother->disease_identity.allele1;
		  
		  r = drand48();
		  if (r < 0.5){
		    allele1 = individual->father->disease_marker.allele2;
		    id1 = individual->father->disease_identity.allele2;
		  }
		  r = drand48();
		  if (r < 0.5){
		    allele2 = individual->mother->disease_marker.allele2;
		    id2 = individual->mother->disease_identity.allele2;
		  }

		  double prob_disease;
		  if (allele1 == minor_allele+1 && allele2 == minor_allele+1)
		    prob_disease = f2;
		  else if (allele1 != minor_allele+1 && allele2 != minor_allele+1)
		    prob_disease = f0;
		  else
		    prob_disease = f1;


		  r = drand48();
		  individual->affection = UNAFFECTED;
		  if (r < prob_disease)
		    {
		      individual->affection = AFFECTED;
		      fam_num_cases++;
		    }

		  individual->disease_marker.allele1 = allele1;
		  individual->disease_marker.allele2 = allele2;
		  individual->disease_identity.allele1 = id1;
		  individual->disease_identity.allele2 = id2;
		  //printf("         NON founder %i] disease alleles (%i,%i) prob_disease %f phenotype %i\n", individual->individual_id, individual->disease_marker.allele1, individual->disease_marker.allele2, prob_disease, individual->affection);	  
		} // end if founder

	    } // end for each individual




	  if (fam_num_cases < required_cases)
	    {
	      isRejected = true;
#ifdef DEBUG
	      printf("         ------------ REJECTED with %i ---------------\n", fam_num_cases);
#endif
	    }
	  else
	    {
	      isRejected = false;
	      number_of_cases += fam_num_cases;
	      printf("         ############ ACCEPTED with %i ###############\n", fam_num_cases);
	    }

	} // end while rejected

    } // end for each family


  printf("Number of cases: %i / %i\n\n", number_of_cases, number_of_individuals);
  return number_of_cases;
}






int sample_genotypes_given_disease_locus(int disease_locus, int minor_allele)
{

  // partition the founder haplotypes into two groups (based on the disease locus)
  // value is the index of a haplotype with the minor allele at the disease locus
  vector<int> minor_haplotypes;
  vector<int> major_haplotypes;

  int number_of_recombinations = 0;

  for (unsigned int i = 0;  i < founder_haplotypes.size();  i++)
    {
      if (founder_haplotypes[i][disease_locus] == minor_allele)
	minor_haplotypes.push_back(i);
      else
	major_haplotypes.push_back(i);
    }
  int count_minor_haps = minor_haplotypes.size();
  int count_major_haps = major_haplotypes.size();


  // For each family
  vector<Family>::iterator fam;
  for (fam = pedigrees.begin();  fam != pedigrees.end();  fam++)
    {

      // for each individual (in top-down order)
      vector<OrderItem>::iterator top_down;
      for (top_down = fam->order.end()-1;  top_down != fam->order.begin()-1;  top_down--)
	{
	  //printf("  (sample phenotypes) Processing indiv %i:%i \n", fam->family_id, top_down->individual_id);
	  //fflush(stdout);
	  
	  
	  Individual* individual = top_down->individual;
	  
	  double r = 0.0;
	  
	  if (individual->father_id == UNKNOWN && individual->mother_id == UNKNOWN)
	    { // is a founder
	      int hap_index1 = -1;
	      int hap_index2 = -1;
	      
	      r = drand48();
	      if (individual->disease_marker.allele1 == minor_allele+1)
		{
		  int hap = (int) floor(r * count_minor_haps);
		  hap_index1 = minor_haplotypes[hap];
		}
	      else
		{
		  int hap = (int) floor(r*count_major_haps);
		  hap_index1 = major_haplotypes[hap];
		}
	      
	      r = drand48();
	      if (individual->disease_marker.allele2 == minor_allele+1)
		{
		  int hap  = (int) floor(r*count_minor_haps);
		  hap_index2 = minor_haplotypes[hap];
		}
	      else
		{
		  int hap = (int) floor(r*count_major_haps);
		  hap_index2 = major_haplotypes[hap];
		}
	      
	      // copy the hapltotypes into the marker vector
	      for (unsigned int a = 0; a < individual->genotypes.size();  a++)
		{
		  individual->genotypes[a].allele1 = founder_haplotypes[hap_index1][a]+1;
		  individual->genotypes[a].allele2 = founder_haplotypes[hap_index2][a]+1;
		  individual->identity[a].allele1 = individual->disease_identity.allele1;
		  individual->identity[a].allele2 = individual->disease_identity.allele2;
		}
	      
	    }
	  else
	    { // is not founder
	      
	      // randomly choose the source (1 or 2) of the maternal allele
	      int mother_haplotype = -1;
	      int m_hap_at_d = -1;
	      Marker dl = individual->disease_marker;
	      //Marker m_disease_locus = individual->mother->disease_marker; 
	      Marker m_disease_locus = individual->mother->genotypes[disease_locus];
	      if (individual->disease_marker.allele2 == m_disease_locus.allele1)
		{
		  if (individual->disease_marker.allele2 == m_disease_locus.allele2)
		    {
		      r = drand48();
		      if (r < 0.5)
			m_hap_at_d = 1;
		      else
			m_hap_at_d = 2;
		    }
		  else
		    m_hap_at_d = 1;
		}
	      else if (individual->disease_marker.allele2 == m_disease_locus.allele2)
		{
		  m_hap_at_d = 2;
		}
	      if (m_hap_at_d == -1)
		{
		  printf("ERROR: child %i:%i (%i,%i) has mendelian inconsistency with mother (%i,%i) at disease locus\n", individual->family_id, individual->individual_id, dl.allele1, dl.allele2, m_disease_locus.allele1, m_disease_locus.allele2);
		  exit(-42);
		}

	      // randomly choose the source (1 or 2) of the paternal allele
	      int father_haplotype = -1;
	      int f_hap_at_d = -1;
	      Marker f_disease_locus = individual->father->genotypes[disease_locus];
	      if (individual->disease_marker.allele1 == f_disease_locus.allele1)
		{
		  if (individual->disease_marker.allele1 == f_disease_locus.allele2)
		    {
		      r = drand48();
		      if (r < 0.5)
			f_hap_at_d = 1;
		      else
			f_hap_at_d = 2;
		    }
		  else
		    f_hap_at_d = 1;
		}
	      else if (individual->disease_marker.allele1 == f_disease_locus.allele2)
		{
		  f_hap_at_d = 2;
		}
	      if (f_hap_at_d == -1)
		{
		  printf("ERROR: child %i:%i (%i,%i) has mendelian inconsistency with father (%i,%i) at disease locus\n", individual->family_id, individual->individual_id, dl.allele1, dl.allele2, f_disease_locus.allele1, f_disease_locus.allele2);
		  exit(-42);
		}


	      // copy the hapltotypes into the marker vector (by two copy operations)
	      // by copying forward from the disease locus 
	      mother_haplotype = m_hap_at_d;
	      father_haplotype = f_hap_at_d;
	      for (unsigned int a = disease_locus; a < individual->genotypes.size();  a++)
		{
		  if (father_haplotype == 1){
		    individual->genotypes[a].allele1 = individual->father->genotypes[a].allele1;
		    individual->identity[a].allele1 = individual->father->identity[a].allele1;
		  } else {
		    individual->genotypes[a].allele1 = individual->father->genotypes[a].allele2;
		    individual->identity[a].allele1 = individual->father->identity[a].allele2;
		  }
		  
		  if (mother_haplotype == 1) {
		    individual->genotypes[a].allele2 = individual->mother->genotypes[a].allele1;
		    individual->identity[a].allele2 = individual->mother->identity[a].allele1;
		  } else {
		    individual->genotypes[a].allele2 = individual->mother->genotypes[a].allele2;
		    individual->identity[a].allele2 = individual->mother->identity[a].allele2;
		  }
		  
		  
		  // recombination ??
		  r = drand48();
		  if (r < recombination_rates[a])
		    {
		      int before = father_haplotype;
		      father_haplotype = 1 - father_haplotype + 2;
		      number_of_recombinations++;
		      printf("    Recombination in paternal haplotype of %i:%i at locus %i from %i to %i\n", individual->family_id, individual->individual_id, a, before, father_haplotype);
		    }
		  r = drand48();
		  if (r < recombination_rates[a])
		    {
		      int before = mother_haplotype;
		      mother_haplotype = 1 - mother_haplotype + 2;
		      number_of_recombinations++;
		      printf("    Recombination in maternal haplotype of %i:%i at locus %i from %i to %i\n", individual->family_id, individual->individual_id, a, before, mother_haplotype);
		    }
		}

	      // and then back from the disease locus
	      mother_haplotype = m_hap_at_d;
	      father_haplotype = f_hap_at_d;
	      for (int a = disease_locus; a >= 0;  a--)
		{
		  if (father_haplotype == 1){
		    individual->genotypes[a].allele1 = individual->father->genotypes[a].allele1;
		    individual->identity[a].allele1 = individual->father->identity[a].allele1;
		  } else {
		    individual->genotypes[a].allele1 = individual->father->genotypes[a].allele2;
		    individual->identity[a].allele1 = individual->father->identity[a].allele2;
		  }
		  
		  if (mother_haplotype == 1){
		    individual->genotypes[a].allele2 = individual->mother->genotypes[a].allele1;
		    individual->identity[a].allele2 = individual->mother->identity[a].allele1;
		  } else {
		    individual->genotypes[a].allele2 = individual->mother->genotypes[a].allele2;
		    individual->identity[a].allele2 = individual->mother->identity[a].allele2;
		  }
		  
		  
		  // recombination ??
		  if (a-1 >= 0)
		    {
		      r = drand48();
		      if (r < recombination_rates[a-1])
			{
			  int before = father_haplotype;
			  father_haplotype = 1 - father_haplotype + 2;
			  number_of_recombinations++;
			  printf("    Recombination in paternal haplotype of %i:%i at locus %i from %i to %i\n", individual->family_id, individual->individual_id, a-1, before, father_haplotype);
			  fflush(stdout);
			}
		      r = drand48();
		      if (r < recombination_rates[a-1])
			{
			  int before = mother_haplotype;
			  mother_haplotype = 1 - mother_haplotype + 2;
			  number_of_recombinations++;
			  printf("    Recombination in maternal haplotype of %i:%i at locus %i from %i to %i\n", individual->family_id, individual->individual_id, a-1, before, mother_haplotype);
			  fflush(stdout);
			}
		    }
		}

	      
	    } // end if founder
	  
	} // end for each individual
    } // end for each family


  printf("\nTotal number of recombinations: %i\n\n", number_of_recombinations);

  return number_of_recombinations;
}




void check_consistency(int disease_locus)
{

  // For each family
  vector<Family>::iterator fam;
  for (fam = pedigrees.begin();  fam != pedigrees.end();  fam++)
    {

      // for non-founder individual (in top-down order)
      vector<OrderItem>::iterator top_down;
      for (top_down = fam->order.end()-1;  top_down != fam->order.begin()-1;  top_down--)
	{
	  //printf("  (sample phenotypes) Processing indiv %i:%i \n", fam->family_id, top_down->individual_id);
	  //fflush(stdout);
	  Individual* individual = top_down->individual;

	  if (individual->genotypes[disease_locus].allele1 != individual->disease_marker.allele1 ||
	      individual->genotypes[disease_locus].allele2 != individual->disease_marker.allele2)
	    {
	      printf("ERROR: inconsistency between sampled genotypes and disease genotype\n");
	      exit(-1);
	    }
	  
	} // end for each individual
    } // end for each family

}



void set_affection(Affection status)
{
    // For each family
  vector<Family>::iterator fam;
  for (fam = pedigrees.begin();  fam != pedigrees.end();  fam++)
    {
      vector<Individual*>::iterator it;
      for (it = fam->members.begin();  it != fam->members.end();  it++)
	{
	  Individual* person = (*it);
	  person->affection = status;
	}
    }
}


int sample_genotypes()
{

  // partition the founder haplotypes into two groups (based on the disease locus)
  // value is the index of a haplotype with the minor allele at the disease locus
  vector<int> minor_haplotypes;
  vector<int> major_haplotypes;

  int number_of_recombinations = 0;

  int count_haps = founder_haplotypes.size();

  int identity_counter = 1;
  
 

  // For each family
  vector<Family>::iterator fam;
  for (fam = pedigrees.begin();  fam != pedigrees.end();  fam++)
    {

      // for non-founder individual (in top-down order)
      vector<OrderItem>::iterator top_down;
      for (top_down = fam->order.end()-1;  top_down != fam->order.begin()-1;  top_down--)
	{
	  //printf("  (sample phenotypes) Processing indiv %i:%i \n", fam->family_id, top_down->individual_id);
	  //fflush(stdout);
	  
	  
	  Individual* individual = top_down->individual;
	  
	  double r = 0.0;
	  
	  if (individual->father_id == UNKNOWN && individual->mother_id == UNKNOWN)
	    { // is a founder
	      r = drand48();
	      int hap_index1 = (int) floor(r * count_haps);
	      
	      r = drand48();
	      int hap_index2  = (int) floor(r*count_haps);
	      
	      // copy the hapltotypes into the marker vector
	      for (unsigned int a = 0; a < individual->genotypes.size();  a++)
		{
		  individual->genotypes[a].allele1 = founder_haplotypes[hap_index1][a]+1;
		  individual->genotypes[a].allele2 = founder_haplotypes[hap_index2][a]+1;
		  individual->identity[a].allele1 = identity_counter;
		  individual->identity[a].allele2 = identity_counter+1;
		}
	      identity_counter+=2;
	      
	    }
	  else
	    { // is not founder
	      
	      int mother_haplotype = -1;
	      r = drand48();
	      if (r < 0.5)
		mother_haplotype = 1;
	      else
		mother_haplotype = 2;

	      int father_haplotype = -1;
	      r = drand48();
	      if (r < 0.5)
		father_haplotype = 1;
	      else
		father_haplotype = 2;

	      // copy the hapltotypes into the marker vector (by two copy operations)
	      // by copying forward from the first locus 
	      for (unsigned int a = 0; a < individual->genotypes.size();  a++)
		{
		  if (father_haplotype == 1){
		    individual->genotypes[a].allele1 = individual->father->genotypes[a].allele1;
		    individual->identity[a].allele1 = individual->father->identity[a].allele1;
		  } else {
		    individual->genotypes[a].allele1 = individual->father->genotypes[a].allele2;
		    individual->identity[a].allele1 = individual->father->identity[a].allele2;
		  }
		  
		  if (mother_haplotype == 1){
		    individual->genotypes[a].allele2 = individual->mother->genotypes[a].allele1;
		    individual->identity[a].allele2 = individual->mother->identity[a].allele1;
		  } else {
		    individual->genotypes[a].allele2 = individual->mother->genotypes[a].allele2;
		    individual->identity[a].allele2 = individual->mother->identity[a].allele2;
		  }
		  
		  
		  // recombination ??
		  r = drand48();
		  if (r < recombination_rates[a])
		    {
		      int before = father_haplotype;
		      father_haplotype = 1 - father_haplotype + 2;
		      number_of_recombinations++;
		      printf("    Recombination in paternal haplotype of %i:%i at locus %i from %i to %i\n", individual->family_id, individual->individual_id, a, before, father_haplotype);
		    }
		  r = drand48();
		  if (r < recombination_rates[a])
		    {
		      int before = mother_haplotype;
		      mother_haplotype = 1 - mother_haplotype + 2;
		      number_of_recombinations++;
		      printf("    Recombination in maternal haplotype of %i:%i at locus %i from %i to %i\n", individual->family_id, individual->individual_id, a, before, mother_haplotype);
		    }
		}

	      
	    } // end if founder
	  
	} // end for each individual
    } // end for each family


  printf("\nTotal number of recombinations: %i\n\n", number_of_recombinations);

  return number_of_recombinations;
}




void print_genotypes(const char* hap_file)
{
  // print out all genotypes (ordered by parental source)
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
	  fprintf(fp, "%i %i ", person->family_id, person->individual_id);
	  fprintf(fp, "%i %i ", person->father_id, person->mother_id);
	  fprintf(fp, "%i %i ", person->sex, person->affection);


	  for (unsigned int snp = 0; snp < person->genotypes.size(); snp++)
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


void print_ibd(const char* hap_file)
{
  // print out all genotypes (ordered by parental source)
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
	  fprintf(fp, "%i %i ", person->family_id, person->individual_id);
	  fprintf(fp, "%i %i ", person->father_id, person->mother_id);
	  fprintf(fp, "%i %i ", person->sex, person->affection);


	  for (unsigned int snp = 0; snp < person->genotypes.size(); snp++)
	    {
	      long allele1 = person->identity[snp].allele1;
	      long allele2 = person->identity[snp].allele2;

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
    }  
}




//////////////////////////////
//
//
int main (int argc, char *argv[])
{
  // read arguements
  if (argc < 11)  // argc is # argvs [1..n] with prog first
    {
      printf ("\nUSAGE: %s [s] [d|n] [prevalence] [risk] [disease allele cutoff] \\\n", argv[0]);
      printf ("               [case freq] [recomb. rates] [hap file] [in pedigree file] \\\n");
      printf ("               [out ped file] (out disease locus) (out ibd file)\n\n");
      printf ("[s]          seed, integer \n");
      printf ("[d|n]        whether to use the disease model, d, or none, n\n");
      printf ("[prevalence] real value, probably somewhere in the range (0.0, 0.4]\n");
      printf ("[risk]       real value >= 1, if risk == 1, then random disease\n");
      printf ("[disease allele cutoff] disease SNP will have minor allele freq. >= cutoff\n");
      printf ("[case frequency]        min freq of affected individuals in the pedigree (rejection sampling)\n");
      printf ("[recombination rates]   file with list of rec. rates separated by white space\n");

      printf ("[hap file]  INPUT list of possible founder haplotypes with alleles {0,1}, format:\n");
      printf ("            snp1 snp2 ... snpN\n");
      printf ("[in pedigree file]   INPUT pedigree (genotypes ignored) in linkage format:\n");
      printf ("            FamID IndID FaID MoID Sex Aff [LOCI ...]\n");
      printf ("[out pedigree file]  OUT pedigree file name\n");
      printf ("(disease locus file) optional OUT file for disease locus which is an integer between [1,num_snps]\n");
      printf ("(ibd file) optional OUT file, but must also have disease file to have ibd output.\n\n");
      printf ("\n");
      printf ("    To simulate the case with no disease locus:\n");
      printf ("         risk = 1\n\n");


      exit(-1);
    }

  const int seed = atof(argv[1]);
  const char use_model = argv[2][0];
  const double prevalence = atof(argv[3]);
  const double risk = atof(argv[4]);
  const double d_freq = atof(argv[5]);
  const double case_freq = atof(argv[6]);
  const char* rec_file = argv[7];
  const char* founder_file = argv[8];
  const char* ped_file = argv[9];
  const char* out_ped_file = argv[10];
  const char* out_dl_file = NULL;
  const char* out_ibd_file = NULL;
  if (argc > 11)
    out_dl_file = argv[11];
  if (argc > 12)
    out_ibd_file = argv[12];


  printf("\nPARAMS: prevalence = %f, risk = %f, d_freq = %f\n", prevalence, risk, d_freq);
  printf("          case freq = %f,\n", case_freq);
  printf("          rec-rates = %s, founders = %s\n", rec_file, founder_file);
  printf("          in ped = %s, out ped = %s, dl file = %s, ibd file = %s\n\n", ped_file, out_ped_file, out_dl_file, out_ibd_file);

  if (d_freq > 0.5)
    {
      printf("FATAL ERROR: There will be no disease allele, since d_freq %f > 0.5\n\n", d_freq);
      exit(-424);
    }



  printf("READING the founder genotypes...\n");
  int number_of_snps = read_founder_file(founder_file);
  read_recombination_file(rec_file, number_of_snps);



  // check feasibility
  int hasReasonableFreq = 0;
  for (int c = 0;  c < number_of_snps;  c++)
    {
      double p = -1.0;
      int minor_allele = -1;
      get_founder_allele_freq(c, p, minor_allele);
      if (p >= d_freq)
	{
	  hasReasonableFreq = 1;
	  c = number_of_snps;
	}
    }
  if (!hasReasonableFreq)
    {
      printf("FATAL ERROR: No founder frequency greater than %f\n", d_freq);
      exit(-66345);
    }


  // seed = time(NULL);
  srand48(seed);

  double p = -1.0; 
  int minor_allele = -1;
  int disease_locus = -1;
  double f0,f1,f2 = -1.1;

  if (use_model == 'd')
    {
      while (p < d_freq)
	{
	  disease_locus = (int) floor(drand48()*number_of_snps);
	  get_founder_allele_freq(disease_locus, p, minor_allele);
	}
      
      
      f0 = prevalence / ((1-p)*(1-p) + 2*p*(1-p)*risk + p*p*risk*risk);
      f1 = risk * f0;
      f2 = risk * risk * f0;
      
      printf("\n");
      printf("disease locus: %i\n", (disease_locus+1));
      printf("disease allele: %i in {1,2},  freq: %f\n", (minor_allele+1), p);
      printf("K: %f, RR: %f, f0 %f, f1 %f, f2 %f\n", prevalence, risk, f0, f1, f2);
      printf("\n");

      // print out all genotypes (ordered by parental source)
      if (out_dl_file != NULL)
      {
        FILE* fp = fopen (out_dl_file, "w");
        if (fp == NULL)
        {	
          printf("Cannot open (w) %s.", out_dl_file);
          exit(-1);
        }
        fprintf(fp, "%i\n", (disease_locus+1));
        fclose(fp);
      }
    }      
      

  printf("READING input pedigree...\n");
  read_pedigree_file(ped_file);
  set_pedigree_pointers();

  printf("set traversal...\n");
  set_traversal_orders();

  printf("set affection...\n");


  if (use_model == 'd')
    {
      // rejection sample the disease
      printf("Sampling phenotypes...\n");
      sample_phenotypes(disease_locus, p, minor_allele, f0, f1, f2, case_freq);
      
      printf("Sampling genotypes conditional on phenotypes...\n");
      allocate_markers(number_of_snps);
      sample_genotypes_given_disease_locus(disease_locus, minor_allele);
      
      check_consistency(disease_locus);
    }
  else
    {
      set_affection(UNAFFECTED);
      allocate_markers(number_of_snps);
      sample_genotypes();
    }


  printf("Printing simulation to file...\n");
  print_genotypes(out_ped_file);

  // Print the IBD information
  if (out_ibd_file != NULL)
    print_ibd(out_ibd_file);


  deleteped();
}


