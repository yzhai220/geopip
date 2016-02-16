
#ifndef IBDSETFOREST
#define IBDSETFOREST

#include <map>
#include <vector>
#include <queue>
#include <set>
#include <iostream>

#include "indicator.h"


#define PLACE_HOLDER -33



#define IBDDEBUG 1
//REMOVE THIS LINE TO ENABLE DEBUG
#undef IBDDEBUG

#ifdef IBDDEBUG
#define debout(x) cout << x << endl
#define deboutnn(x) cout << x << flush
#define deboutln(x) cout << x << endl << flush
#define deb(x) x
#else
#define debout(x) ;
#define deboutln(x) ; 
#define deboutnn(x) ; 
#define deb(x)  ;
#endif




/////////////////
// Set forest for integers numbered [1,2,...,size]
//
// Invarient condition:
//   flat_set_forest[i] <= i
//   tall_set_forest[i] <= i
//
// flat_set_forest[i] and tall_set_forest[i] have same root element
//

class IBDSetForest
{
public:
  IBDSetForest(){ 
    total_children = 0;
    setalleles_result = true;
    num_cc = 0;
  }  

  // number of people === _size
  IBDSetForest(int _size, int _num_cc){
    flat_set_forest = vector<int>(2*_size, -1);
    tall_set_forest = vector<int>(2*_size, -1);

    genotypes = vector<int>(2*_size, 0);
    set_alleles = vector<int>(4*_size, 0);
    num_children = vector<int>(_size, 0);



    alleles1 = vector<int>(2*_size, -1);
    alleles2 = vector<int>(2*_size, -1);

    which_cc = vector<int>(2*_size,-1);

    num_cc = _num_cc;


    total_children = 0;
    setalleles_result = true;
    num_alleles = 0;
  }  
  IBDSetForest (const IBDSetForest& c);


  // allele1 & allele2 >= 0, where 0 indicates ungenotyped
  int InsertNode(int node, int ped_index, int allele1, int allele2, int _num_children);

  void SetNumAlleles(int _num_alleles){ num_alleles = _num_alleles; }



  // Find the set this node (2*(_node-1)+allele) belongs to
  // returns -1 for error or [0..size-1] for the set
  int FindSet(int _node, int allele);


  // takes pedigree index for individual  [0..n]
  int GetAlleleIndicatorPed(int index, int parent); // index starting from zero
  int GetAlleleIndicator(int indiv, int allele, int parent); // indiv in our order starting from one

  // Return the list of alleles for a particular set
  vector<int> GetAlleles(int indiv, int allele);

  // return a weight for the inheritance path based on given weights
  double GetLogEmissionWeight(vector<Indicator>& indicator_weights);

  // return a probability for observed genotypes given this inheritance path
  double GetLogEmissionProb();

  // return a list of haplotype orientations 
  // for the IBD inheritance path that we have here
  vector<int> GetOrientations();

  // return a list of alleles for every person where
  // the founder alleles are sampled uniformly from the 1-2 options per CC
  vector<int> SampleAlleles();
  vector<int> SampleInfiniteAlleles();
  double GetLogSampleProb(){return log_sample_prob;};


  // Find the next node having the same greatest ancestor (belonging to same set)
  // prev element starts at ZERO and subsequently is previous return-value of NextSetElement
  // returns -1 for error, 0 for finished, node number otherwise
  int NextSetElement(int ancestor, int prev_element);

  // Find the next node having the same clade ancestor (belonging to same sub-set)
  // prev element starts at ZERO and subsequently is previous 
  // return-value of NextSubSetElement
  // returns -1 for error, 0 for finished, node number otherwise
  int NextSubSetElement(int ancestor, int prev_element);




  // Perform Union on the two sets
  // Post-conditions: 
  //   1) flat_set_forest is flat
  //   2) flat_set_forest[i] >= i
  // Simultaneously perform union with both alleles of u to the
  // specified alleles of parents v and w.
  // Checks for consistency and returns either true or false.
  bool Union(int u, int v, bool allele_v, int w, bool allele_w);


  void UnionWithoutSetAlleles(int u, int v, bool allele_v, int w, bool allele_w);

  bool SetAlleles();


  // Pre-condition:
  //   1) flat_set_forest[i] >= i
  //   2) tall_set_forest[i] contains memory of union history (i.e. the subsets)
  // Post-condition:
  //   1) flat_set_forest and tall_set_forest are returned to the state they had before
  //      element u joined any of it's ancestor sets
  //   2) flat_set_forset[i] >= i
  int ReverseUnion(int u);


  void PrintVectors();


  // Perform Union on the two sets
  // Post-conditions: 
  //   1) flat_set_forest is flat
  //   2) flat_set_forest[i] >= i
  // returns bool indicating whether the inheritance path has valid allele settings
  int Union(int u, bool allele_u, int v, bool allele_v);


  bool validAlleleSets(){return setalleles_result;};


protected:

  // Find the set this node belongs to
  // returns -1 for error or [0..size-1] for the set
  int FindSet(int node);

  bool IsRoot(int node);



  unsigned int GetIndiv(int index){
    return (unsigned int) floor((double) index / 2);
  }


  void CleanUp(){
    for (unsigned int i = 0;  i < set_alleles.size();  i++)
      set_alleles[i] = 0;
    for (unsigned int i = 0;  i < alleles1.size();  i++)
      alleles1[i] = -1;
    for (unsigned int i = 0;  i < alleles2.size();  i++)
      alleles2[i] = -1;
    for (unsigned int i = 0;  i < which_cc.size();  i++)
      which_cc[i] = -1;

  }
 

  bool recursiveSetAlleles(unsigned int index, 
			   int allele1, 
			   vector<int>& alleles,
			   vector<int>& visited);
  void markCC(unsigned int index, int cc, unsigned int& genotyped);
  void maskCC(int cc, vector<int>& alleles, int value);
    
  // Return the list of alleles for a particular set
  vector<int> GetAlleles(int index);


  vector<int> flat_set_forest;
  vector<int> tall_set_forest;

  vector<int> num_children;
  vector<int> genotypes;
  vector<int> set_alleles;


  vector<int> alleles1;
  vector<int> alleles2;


  // == 2*(# non-founders), since each child is counted 
  // twice, once by mother, once by father
  int total_children; 

  bool setalleles_result;

  int num_cc;
  int num_alleles;  // use to compute prob. of observation given inheritance path

  // after SetAlleles()
  //   cc_alleles has 2 alleles per person (cc_allele1 & cc_allele2)
  //vector<int> cc_alleles;
  vector<int> which_cc;

  double log_sample_prob;

  // maping from data structure index to pedigree index
  map<int, int> indices;
};








IBDSetForest::IBDSetForest (const IBDSetForest& c)
{
  indices = c.indices;
  if (c.flat_set_forest.size() != flat_set_forest.size())
    {
      printf("WARNING: size mismatch in IBDSetForest::IBDSetForest\n");
    }

  for (unsigned int i = 0; i < c.flat_set_forest.size();  i++)
    {
      flat_set_forest[i] = c.flat_set_forest[i];
    }
  for (unsigned int i = 0; i < c.tall_set_forest.size();  i++)
    {
      tall_set_forest[i] = c.tall_set_forest[i];
    }

  for (unsigned int i = 0; i < c.genotypes.size();  i++)
    {
      genotypes[i] = c.genotypes[i];
    }
  for (unsigned int i = 0; i < c.set_alleles.size();  i++)
    {
      set_alleles[i] = c.set_alleles[i];
    }
  for (unsigned int i = 0; i < c.num_children.size();  i++)
    {
      num_children[i] = c.num_children[i];
    }

  for (unsigned int i = 0; i < c.alleles1.size();  i++)
    {
      alleles1[i] = c.alleles1[i];
    }
  for (unsigned int i = 0; i < c.alleles2.size();  i++)
    {
      alleles2[i] = c.alleles2[i];
    }


  for (unsigned int i = 0; i < c.which_cc.size();  i++)
    {
      which_cc[i] = c.which_cc[i];
    }


  num_cc = c.num_cc;

  total_children = c.total_children;
  setalleles_result = c.setalleles_result;
  num_alleles = c.num_alleles;
  
}









bool IBDSetForest::SetAlleles()
{
  setalleles_result = false;
  CleanUp();

  int cc = 0;
  bool unprocessed = 1;
  while(unprocessed)
    {
      deboutln("CC: "<< cc << " unprocessed: " << unprocessed);

      // find the first unprocessed node (i.e. allele in flat_set_forest)
      unsigned int idx = 0; // index in flat_set_forest
      while (which_cc[idx] != -1)
	{
	  idx++;
	  if (idx >= which_cc.size())
	    break;
	}
      if (idx >= which_cc.size())
	break;
      deboutln("  found index "<<idx);


      // also finds a genotyped individual in the CC
      markCC(idx, cc, idx);


      deboutln("     mark found geno "<<idx<<" for cc "<< cc);


      // do recursion 
      // only if we have a genotyped person
      int indiv = GetIndiv(idx);
      int a1 = genotypes[2*indiv];
      int a2 = genotypes[2*indiv+1];
      //printf("  running on typed indiv: %i, with geno %i,%i\n", idx, a1, a2);

      //if (a1 != 0 || a2 != 0)
      //{
	  vector<int> visit1 = vector<int>(alleles1.size(), 0);
	  bool compat_a1 = recursiveSetAlleles(idx, a1, alleles1, visit1);
	  vector<int> visit2 = vector<int>(alleles1.size(), 0);
	  bool compat_a2 = recursiveSetAlleles(idx, a2, alleles2, visit2);

	  if (!compat_a1 && !compat_a2)
	    {
	      CleanUp();
	      return false;
	    }

	  if (!compat_a1)
	    maskCC(cc, alleles1, -1);
	  if (!compat_a2)
	    maskCC(cc, alleles2, -1);
	  //}

      cc++;
    } // while unprocessed
  num_cc = cc;

  // Copy results over to set_alleles
  for (unsigned int i = 0;  i < alleles1.size();  i++)
    {
      int indiv = GetIndiv(i);
      if (alleles1[i] == 0 && genotypes[2*indiv] != 0)
	printf("WARNING: index %i (indiv %i) is genotyped but has all alleles\n", i, indiv);
      set_alleles[2*i] = alleles1[i];
    }
  for (unsigned int i = 0;  i < alleles2.size();  i++)
    {
      int indiv = GetIndiv(i);
      if (alleles2[i] == 0 && genotypes[2*indiv] != 0)
	printf("WARNING: index %i (indiv %i) is genotyped but has all alleles\n", i, indiv);
      set_alleles[2*i+1] = alleles2[i];
    }

  //PrintVectors();

  setalleles_result = true;
  return true;
}


void IBDSetForest::maskCC(int cc, vector<int>& alleles, int value)
{
  for (unsigned int i = 0;  i < which_cc.size();  i++)
    if (which_cc[i] == cc)
      alleles[i] = value;
}


// Mark all the individuals in this connected component
// set genotyped to the index of some genotyped allele
void IBDSetForest::markCC(unsigned int idx, int cc, unsigned int& genotyped)
{
  if (which_cc[idx] == -1)
    {
      //deboutln("     mark which_cc ["<<idx<<"] = "<< cc);
      which_cc[idx] = cc;

      // geno compatible edges
      unsigned int indiv = GetIndiv(idx);
      if (genotypes[2*indiv] != 0 && genotypes[2*indiv+1] != 0)
	{
	  if (2*indiv != idx)
	    markCC(2*indiv, cc, genotyped);
	  if (2*indiv+1 != idx)
	    markCC(2*indiv+1, cc, genotyped);
	}

 
      // Identity edges (all other alleles in this set)
      unsigned int ancestor = FindSet(idx);
      int prev_element = -1;
      int element = NextSetElement(ancestor, prev_element);
      while (element > -1)
	{
	  if ((unsigned) element != idx)
	    {
	      //deboutln("        element "<< element);

	      markCC(element, cc, genotyped);
	    }
	  prev_element = element;
	  element = NextSetElement(ancestor, prev_element);
	} // end all elements

      // set genotyped index to first encountered genotyped person
      if (genotypes[2*indiv] != 0 && genotypes[2*indiv+1] != 0)
	genotyped = 2*indiv;

    }

}



bool IBDSetForest::recursiveSetAlleles(unsigned int idx, 
				       int a1, 
				       vector<int>& alleles,
				       vector<int>& visited)

{
  deboutln("  recursiveSetAlleles("<<idx<<", "<<a1<<", ...)");

  if (alleles[idx] == -1){
    alleles[idx] = a1;
  } else {
    if (alleles[idx] == 0 && a1 != 0)
      printf("WARNING: allele %i was assigned any allele before a specific allele\n", idx);
    if (alleles[idx] != a1)
      return false;
  }

  if (visited[idx] == 1)
    {
      // check genotype
      //int indiv = GetIndiv(idx);
      //if (alleles[idx] == 0 && genotypes[2*indiv] != 0)
      //printf("WARNING: index %i (indiv %i) is genotyped but has all alleles\n", idx, indiv);
      return true;
    }
  visited[idx] = 1;


  // Genotype compatible edges
  unsigned int indiv = GetIndiv(idx);

  unsigned int other_index = 2*indiv;
  if (2*indiv == idx)
    other_index = 2*indiv+1;

  int other_a = -1;
  int b1 = genotypes[2*indiv];
  int b2 = genotypes[2*indiv+1];

  //if (b1 == 0 || b2 == 0)
  //return true;

  if (b1 != 0 && b2 != 0) // otherwise skip to identity edges
    {
      if (a1 == b1)
	other_a = b2;
      if (a1 == b2)
	other_a = b1;
      if (other_a == -1){
	alleles[idx] = -1;
	return false;
      }
      
      bool success = recursiveSetAlleles(other_index, other_a, alleles, visited);
      if (!success)
	return false;
    } 


  // Identity edges (all other alleles in this set)
  unsigned int ancestor = FindSet(idx);
  int prev_element = -1;
  int element = NextSetElement(ancestor, prev_element);
  while (element > -1)
    {
      if ((unsigned) element != idx)
	{
	  deboutln("        element "<< element);

	  bool success = recursiveSetAlleles(element, a1, alleles, visited);
	  if (!success)
	    return false;

	}
      prev_element = element;
      element = NextSetElement(ancestor, prev_element);
    } // end all elements

  return true;
} // end recursiveSetAlleles












// allele1 & allele2 >= 0, where 0 indicates ungenotyped
int IBDSetForest::InsertNode(int node, int ped_index, int allele1, int allele2, int _num_children){
    node--;
    if (node < 0 || 2*node+1 > (signed) flat_set_forest.size()-1)
      {
	printf("ERROR: node %i is out of bounds\n", node);
	exit(-1);
      }
    if (allele1 < 0 || allele2 < 0)
      {
	printf("ERROR: invalid genotype (%i,%i), must be at least zero\n", allele1,allele2);
	exit(-1);
      }

    genotypes[2*node] = allele1;
    genotypes[2*node+1] = allele2;

    num_children[node] = _num_children;
    
    total_children += _num_children;

    // set mapping
    indices[ped_index] = node;


    return 1; // added
}




// Find the set this node belongs to
// returns -1 for error or [0..size-1] for the set
int IBDSetForest::FindSet(int node){
    if (node < 0 || node > (signed) flat_set_forest.size()-1){
      printf("ERROR: FindSet(%i) out of bounds ancestor\n", node);
      return -1;
    }
    // find the geatest ancestor of a node
    int ancestor = node;
    // root has -1, PLACE_HOLDER, or self-pointer
    while (flat_set_forest[ancestor] != -1 && 
	   flat_set_forest[ancestor] != PLACE_HOLDER &&
	   flat_set_forest[ancestor] != ancestor){  
      ancestor = flat_set_forest[ancestor];  // pointers are to parents
    }
    flat_set_forest[ancestor] = ancestor;
    return ancestor;
}

// returns true if the node is a founder/root
bool IBDSetForest::IsRoot(int node)
{
  if (node == FindSet(node))
    return true;
  return false;
}




// Find the set this node belongs to
// returns -1 for error or [0..size-1] for the set
int IBDSetForest::FindSet(int _node, int allele){
  int node = 2*(_node-1)+allele;

  if (node < 0 || node > (signed) flat_set_forest.size()-1){
    printf("ERROR: FindSet(%i,%i) out of bounds ancestor\n", _node,allele);
    return -1;
  }
  // find the geatest ancestor of a node
  int ancestor = node;
  // root has -1, PLACE_HOLDER, or self-pointer
  while (flat_set_forest[ancestor] != -1 && 
	 flat_set_forest[ancestor] != PLACE_HOLDER &&
	 flat_set_forest[ancestor] != ancestor){  
    ancestor = flat_set_forest[ancestor];  // pointers are to parents
  }
  flat_set_forest[ancestor] = ancestor;
  return ancestor;
}


///////
// iterate through indicator_weights and sum the log contributions of
// the weights appearing in this inheritance path
double IBDSetForest::GetLogEmissionWeight(vector<Indicator>& indicator_weights)
{
  double weight = 0;

  //printf("indicator.size(): %i\n", (int) indicator_weights.size());
  //fflush(stdout);

  for (unsigned int i = 0;  i < indicator_weights.size();  i++)
    {
      Indicator *indi = &indicator_weights[i];
      int index = indi->indiv_index;
      int parent = indi->parent;

      int node = indices[index]; // node in the inheritance path
      int inheritance_pointer = flat_set_forest[2*node+parent];
      if (inheritance_pointer % 2 == 0)
  	weight += log(indi->grandfather);
      else
  	weight += log(indi->grandmother);
    }

  return weight;
}


/////////////
//
// Return probability for observed genotypes given this inheritance path
// Use cc_alleles and num_alleles to compute uniform probabilties.
// 
double IBDSetForest::GetLogEmissionProb()
{
  double unif_allele_freq = 1.0 / (double) num_alleles;

  double log_prob = 0;

  // independent probabilities for each connected component
  // there are two disjoint probabilities for each CC (alleles1, alleles2)
  vector<double> cc_prob0(num_cc, 1.0);
  vector<double> cc_prob1(num_cc, 1.0);

  //PrintVectors();


  deboutnn("GetLogEmissionProb: ");    

  // take product of cc founder contributions
  // for cc_prob0 and cc_prob1
  vector<bool> isDifferent(num_cc, 0);
  for (unsigned int i = 0;  i < flat_set_forest.size();  i++)
    {
      // if we have a root
      if (flat_set_forest[i] == (signed) i || flat_set_forest[i] == -1){
	vector<int> alleles = GetAlleles(i);
	deboutnn(alleles[0] << "," << alleles[1]<< "; ");

	// Pr[allele] = {0, p_a, 1}
	//   if no allele, allele a, or any allele, respectively

	if (alleles[0] == 0 || alleles[1] == 0)
	  {
	    continue;
	  }
	if (alleles[0] != alleles[1])
	  {
	    isDifferent[which_cc[i]] = 1;
	  }

	assert (alleles[0] != -1 || alleles[1] != -1);

	//if (alleles[0] != 1)
	  cc_prob0[which_cc[i]] *= unif_allele_freq;
	//if (alleles[1] != 1)
	  cc_prob1[which_cc[i]] *= unif_allele_freq;
      }
    }
  deboutln("");


  // product over cliques of cc_prob0[c] + cc_prob1[c]
  for (unsigned int c = 0;  c < cc_prob0.size();  c++)
    {
      /* double sum = 0.0; */
      /* if (cc_prob0[c] > 0.0 && cc_prob1[c] > 0.0) */
      /* 	sum = 0.5 * cc_prob0[c] + 0.5*cc_prob1[c]; */
      /* else */
      /* 	sum = cc_prob0[c] + cc_prob1[c]; */

      double sum = 0.0;
      if (!isDifferent[c])
	sum = cc_prob0[c];
      else
	sum = cc_prob0[c] + cc_prob1[c];


      //printf("sum: %f \n", sum);

      if (sum < 0 || 1 < sum)
	printf("WARNING: emission probability is not a probability\n");
      log_prob += log(sum);
    }


  //if (log_prob <= -999999)
  //{
      deboutln("prob: " << exp(log_prob));
      //}

  return log_prob;
}


// parent 0 or 1
// takes pedigree index for individual  [0..n]
// ped index is order individuals appear in ped file
int IBDSetForest::GetAlleleIndicatorPed(int index, int parent)
{
  int node = indices[index];
  return flat_set_forest[2*node+parent] % 2;
}

// take indiv number in data structure order
int IBDSetForest::GetAlleleIndicator(int indiv, int allele, int parent)
{
  int ancestor = 2*(indiv-1) + allele;
  return flat_set_forest[2*ancestor+parent] % 2;
}


vector<int> IBDSetForest::GetAlleles(int indiv, int allele){
  int ancestor = 2*(indiv-1) + allele;
  return GetAlleles(ancestor);
}


// Return the list of possible alleles for a particular individual and allele
//
// In order to get only possible alleles, need to use cc_alleles
// and mask over the set_alleles for this person (to mask out impossible alleles).
//
vector<int> IBDSetForest::GetAlleles(int index){
  vector<int> rtn = vector<int>(2, -1);

  if (index < 0 || index > (signed) flat_set_forest.size()-1){
    printf("ERROR: GetAlleles(%i) out of bounds indiv\n", index);
    return rtn;
  }


  if (set_alleles[2*index] != -1)
    {
      rtn[0] = set_alleles[2*index];
    }
  if (set_alleles[2*index+1] != -1)
    {
      rtn[1] = set_alleles[2*index+1];
    }
  return rtn;
}



// Find the next node having the same greatest ancestor (belonging to same set)
// prev element starts at ZERO and subsequently is previous return-value of NextSetElement
// returns -1 for error, 0 for finished, node number otherwise
int IBDSetForest::NextSetElement(int ancestor, int prev_element){
  if (ancestor < 0 || ancestor > (signed)flat_set_forest.size()-1){
    printf("ERROR: NextSetElement(%i, %i) out of bounds ancestor\n", ancestor, prev_element);
    return -1;
  }
  if (prev_element >= (signed)flat_set_forest.size()-1){
    return -1; // finished
  }
  int this_element = prev_element; // an already checked location
  int this_ancestor = -1;
  while (this_ancestor != ancestor){
    this_element++;
    if (this_element > (signed)flat_set_forest.size()-1){
      return -1;
    }
    this_ancestor = FindSet(this_element);
    //printf("  FindSet(%i) = %i  cmp w/ %i\n", this_element, this_ancestor, ancestor);
  }
  return this_element;
}








// Find the next node having the same clade ancestor (belonging to same sub-set)
// prev element starts at ZERO and subsequently is previous return-value of NextSetElement
// returns -1 for error, 0 for finished, node number otherwise
int IBDSetForest::NextSubSetElement(int ancestor, int prev_element){
  if (ancestor < 0 || ancestor > (signed)tall_set_forest.size()-1){
    printf("ERROR: NextSubSetElement(%i, %i) out of bounds ancestor\n", ancestor, prev_element);
    return -1;
  }
  if (prev_element >= (signed)tall_set_forest.size()-1){
    return -1; // finished
  }
  int this_element = prev_element; // an already checked location
  int this_ancestor = -1;
  while (this_ancestor != ancestor){
    this_element++;
    if (this_element > (signed)tall_set_forest.size()-1){
      return -1;
    }
    this_ancestor = tall_set_forest[this_element];
  }
  return this_element;
}







// Perform Union on the two sets
// Post-conditions: 
//   1) flat_set_forest is flat
//   2) flat_set_forest[i] >= i
int IBDSetForest::Union(int u, bool allele_u, int v, bool allele_v){
  int ancestor_u = FindSet(u, allele_u);
  int ancestor_v = FindSet(v, allele_v);

  //printf("Union(u=%i,au=%i, v=%i,av=%i)\n",u,allele_u,v,allele_v);
  //printf(" ancestor_u %i, ancestor_v %i\n",ancestor_u, ancestor_v);

  u--;
  v--;
  int node_u = 2*u+allele_u;
  int node_v = 2*v+allele_v;

  // already joined
  if (ancestor_u == ancestor_v)
    return 2;


  // condition maintains the invarient tall_set_forest[i] >= i
  if (u <= v)
    tall_set_forest[node_u] = node_v;
  else 
    tall_set_forest[node_v] = node_u;



  // condition maintains the invarient flat_set_forest[i] >= i
  if (ancestor_u > ancestor_v){ 
    // change all elements of Set(v) to point to root of u
    int prev_element = -1;
    int next_element = NextSetElement(ancestor_v, prev_element);
    while (next_element > -1){
      flat_set_forest[next_element] = ancestor_u;
      prev_element = next_element;
      next_element = NextSetElement(ancestor_v, prev_element);
    }
    flat_set_forest[ancestor_v] = ancestor_u;

  } else { // ancestor_v > ancestor_u
    // change all elements of Set(u) to point to root of v
    int prev_element = -1;
    int next_element = NextSetElement(ancestor_u, prev_element);
    while (next_element > -1){
      flat_set_forest[next_element] = ancestor_v;
      prev_element = next_element;
      next_element = NextSetElement(ancestor_u, prev_element);
    }
    flat_set_forest[ancestor_u] = ancestor_v;
  }

  return 1;
}




// Perform Union on the two sets
// Post-conditions: 
//   1) flat_set_forest is flat
//   2) flat_set_forest[i] >= i
// Simultaneously perform union with both alleles of u to the
// specified alleles of parents v and w.
// Checks for consistency and returns either true or false.
bool IBDSetForest::Union(int u, int v, bool allele_v, int w, bool allele_w){

  Union(u,0,v,allele_v);
  Union(u,1,w,allele_w);




  u--;
  v--;
  w--;




  CleanUp();
  bool rtn = SetAlleles();

  //if (!rtn)
  //printf(" ####### failed #############\n");

  return rtn;

}






// Perform Union on the two sets
// Post-conditions: 
//   1) flat_set_forest is flat
//   2) flat_set_forest[i] >= i
// Simultaneously perform union with both alleles of u to the
// specified alleles of parents v and w.
// Checks for consistency and returns either true or false.
void IBDSetForest::UnionWithoutSetAlleles(int u, int v, bool allele_v, int w, bool allele_w){

  Union(u,0,v,allele_v);
  Union(u,1,w,allele_w);

  u--;
  v--;
  w--;
}








// Pre-condition:
//   1) flat_set_forest[i] >= i
//   2) tall_set_forest[i] contains memory of union history (i.e. the subsets)
// Post-condition:
//   1) flat_set_forest and tall_set_forest are returned to the state they had before
//      element u joined any of it's ancestor sets
//   2) flat_set_forset[i] >= i
int IBDSetForest::ReverseUnion(int u)
{
  u--;
  int node_u = 2*u;
  int node_u2 = 2*u+1;


  //int ancestor_u = FindSet(node_u);

  if (flat_set_forest[node_u] == node_u && flat_set_forest[node_u2] == node_u2)
    return 2;

  tall_set_forest[node_u] = node_u;
  flat_set_forest[node_u] = node_u;
  tall_set_forest[node_u2] = node_u2;
  flat_set_forest[node_u2] = node_u2;

  // u,0
  int prev_element = -1;
  int next_element = NextSubSetElement(node_u, prev_element);
  while (next_element > -1){
    flat_set_forest[next_element] = node_u;
    prev_element = next_element;
    next_element = NextSubSetElement(node_u, prev_element);
  }

  // u,1
  prev_element = -1;
  next_element = NextSubSetElement(node_u2, prev_element);
  while (next_element > -1){
    flat_set_forest[next_element] = node_u2;
    prev_element = next_element;
    next_element = NextSubSetElement(node_u2, prev_element);
  }



  CleanUp();

  //return SetAlleles();
  return 1;
}









void IBDSetForest::PrintVectors()
{
  printf("------------------\n");

  printf("Flat Set Forest: \n");
  for(unsigned int u = 0; u < flat_set_forest.size()-1;  u+=2)
    printf(" %i] %i,%i ", u/2+1, flat_set_forest[u], flat_set_forest[u+1]);
  //printf(" %i] %i,%i|%i,%i ", u/2+1, flat_set_forest[u]/2, flat_set_forest[u]%2, flat_set_forest[u+1]/2,flat_set_forest[u+1]%2);

  printf("\nTall Set Forest: \n");
  for(unsigned int u = 0; u < tall_set_forest.size()-1;  u+=2)
    printf(" %i] %i,%i ", u/2+1, tall_set_forest[u], tall_set_forest[u+1]);

  printf("\n\nTotal Children: %i", total_children);
  printf("\nNum Children: \n");
  for(unsigned int u = 0; u < num_children.size();  u++)
    printf(" %i] %i ", u+1, num_children[u]);

  printf("\nGenotypes: \n");
  for(unsigned int u = 0; u < genotypes.size()-1;  u+=2)
    printf(" %i] %i,%i ", u/2+1, genotypes[u], genotypes[u+1]);
      
  printf("\n\nSet Alleles: \n");
  for(unsigned int u = 0; u < set_alleles.size()-1;  u+=2)
    printf(" %i] %i,%i ", u/2+1, set_alleles[u], set_alleles[u+1]);
 
  printf("\n\nWhich CC: \n");
  for (unsigned int i = 0;  i < which_cc.size();  i+=2)
    printf(" %i] %i,%i ", i/2+1, which_cc[i], which_cc[i+1]);


  printf("\n------------------\n");
}




vector<int> IBDSetForest::GetOrientations()
{
  // Orientations only on typed individuals
  // consider in the order they are in this set-forest

  // MUST consider the connected components, because
  // the orientations (i.e. paternal-maternal origin)
  // of the haplotypes in those individuals are
  // tied together.

  // Within a connected component, we don't have symmetry.
  // Between connected components, we have symmetries.

  

  //PrintVectors();

  vector<vector<int> > cc_orient0(num_cc, vector<int>(0,0));
  vector<vector<int> > cc_orient1(num_cc, vector<int>(0,0));


  // find the orientation of all the genotyped individuals
  int base = 1;
  for (unsigned int i = 0;  i < flat_set_forest.size();  i+=2, base*=2)
    {

      // recall that geno phasings are in columns here (i.e. alleles0[0] & alleles1[0])
      vector<int> alleles0 = GetAlleles(i);
      vector<int> alleles1 = GetAlleles(i+1);

      
      // given genotype
      int g1 = genotypes[i];
      int g2 = genotypes[i+1];
      // ignore untyped individuals
      if (g1==0 && g2==0)
	continue;

      // i and it's neighbor allele have same cc
      // since i is genotyped
      int cc = which_cc[i]; 
      
      //printf("alleles %i,%i cc %i with geno (%i,%i) and alleles constraints (%i,%i,) or (%i,%i)\n", i, i+1, cc, g1,g2, alleles0[0],alleles1[0], alleles0[1],alleles1[1]);

      
      // Alleles can contain numbered allele > 0
      //    or 0 for any allele or -1 for no allele
      // Since these alleles must explain the genotype,
      //  - if alleles0[i] == -1, then alleles1[i] == -1
      //  - if alleles0[i] == g1, then alleles2[i] == g2
      //  - if alleles0[i] == g2, then alleles2[i] == g1
      //  - none of them should be zero
      if ((alleles0[0] == -1 && alleles1[0] != -1) || 
	  (alleles0[1] == -1 && alleles1[1] != -1))
	{
	  printf("Error in masking of column: (%i,%i,) or (%i,%i)\n",alleles0[0],alleles1[0], alleles0[1],alleles1[1]);
	  exit(-5967643);
	}
      if (alleles0[0] == 0 || alleles0[1] == 0 || alleles1[0] == 0 || alleles1[1] == 0)
	{
	  printf("Error: there should be no zero in the GetAllele result for a typed individual\n");
	  exit(-5967644);
	}
      
      // Compute orientation as:
      // cc_orient#[cc] += I{oriented opposite to geno-allele order} * base
      // So, same order has zero contribution

      // if genotype is homozygous, then we need to branch and consider both orientations
      // (in both cc_orient0 and cc_orient1

      int orient0 = -1;      
      int orient1 = -1;      
      if (alleles0[0] == g2 && alleles1[0] == g1)
	orient0 = base;  // oriented opposite to geno-allele order
      if (alleles0[1] == g2 && alleles1[1] == g1)
	orient1 = base;  // oriented opposite to geno-allele order

      if (orient0 != -1)
	{
	  if (cc_orient0[cc].size() == 0)
	    cc_orient0[cc].push_back(0);
	  int size0 = cc_orient0[cc].size();
	  // process only indexes that are already present
	  for (int j = 0;  j < size0; j++)
	    {
	      int index = cc_orient0[cc][j];
	      cc_orient0[cc][j] = index + orient0;
	      if (g1 == g2)	  // branch if homozygous
		cc_orient0[cc].push_back(index + base - orient0);
	    }
	}

      if (orient1 != -1)
	{
	  if (cc_orient1[cc].size() == 0)
	    cc_orient1[cc].push_back(0);
	  int size1 = cc_orient1[cc].size();
	  // process only indexes that are already present
	  for (int j = 0;  j < size1; j++)
	    {
	      int index = cc_orient1[cc][j];
	      cc_orient1[cc][j] = index + orient1;
	      if (g1 == g2)	  // branch if homozygous
		cc_orient1[cc].push_back(index + base - orient1);
	    }
	}

    } // end for each typed person


  // Make one list, and remove duplicates
  // i.e. elements of cc_orient1[cc] that is cc_orient0[cc] 
  // Have one set for each cc:
  vector<set<int> > cc_orient = vector<set<int> >(num_cc, set<int>());
   for (unsigned int cc = 0;  cc < cc_orient0.size();  cc++)
    {
      for (unsigned int j = 0;  j < cc_orient0[cc].size(); j++)
	cc_orient[cc].insert(cc_orient0[cc][j]);
      for (unsigned int j = 0;  j < cc_orient1[cc].size(); j++)
	cc_orient[cc].insert(cc_orient1[cc][j]);
    }



  // Consider the symmetries between the CC's
  // i.e. for all combinations of {0,1} orientations for every CC,
  //      sum the cc orientation integers
  
  vector<int> partial_indices;
  partial_indices.push_back(0);
  for (unsigned int cc = 0;  cc < cc_orient0.size();  cc++)
    {
      int pi_size = partial_indices.size();
      // process only indexes that are already present
      for (int j = 0;  j < pi_size; j++)
	{
	  // we need to branch: seperately sum each option
	  // into the partial orientation index

	  int index = partial_indices[j];
	  if (!cc_orient[cc].empty())
	    {
	      set<int>::iterator si = cc_orient[cc].begin();
	      partial_indices[j] = index + *si;
	      for (si++; si != cc_orient[cc].end(); si++)
		partial_indices.push_back(index + *si);
	    }

	}
    } // end for each cc, enumeration loop


  // check for duplicate orientation indexes


  // Now, partial_indicies list contains all the (full) indices we are interested in.
/*   printf("Orientations: "); */
/*   for (unsigned int i = 0; i < partial_indices.size();  i++) */
/*     printf(" %i", partial_indices[i]); */
/*   printf("\n"); */

  return partial_indices;
}




// pick CC alleles uniformly from the 1-2 options
// if any allele possible, give it the major allele (i.e. 1)
vector<int> IBDSetForest::SampleAlleles()
{
  log_sample_prob = 0.0;


  // for each CC
  vector<int> sampled_cc(num_cc, -1);
  for (unsigned int i = 0;  i < which_cc.size();  i++)
    {
      int cc = which_cc[i];
      if(sampled_cc[cc] == -1)
	{
	  vector<int> alleles = GetAlleles(i);
	  //printf("cc %i] sampling for allele %i from %i,%i... ", cc, i+1, alleles[0], alleles[1]);
	  // if any allele
	  if (alleles[0] == 0 || alleles[1] == 0)
	    {
	      sampled_cc[cc] = 0;
	    }
	  else // otherwise
	    {

	      int num_alleles = 0;
	      for (unsigned int a = 0;  a < 2;  a++){
		if (alleles[a] != -1){
		  num_alleles++;
		}
	      }

	      //printf("num alleles %i...", num_alleles);

	      double u = drand48();
	      double value = 1/(double) num_alleles;
	      log_sample_prob += log(value);
	      for (int a = 0;  a < 2; a++)
		{
		  if (alleles[a] != -1){
		    if (u < value)
		      {
			sampled_cc[cc] = a;
			break;
		      }
		  }
		  value = 1;
		}

	    }
	  //printf(" ...sampled orientation %i\n", sampled_cc[cc]);
		
	}
    }

  // copy selected alleles into an allele vector, one entry per pedigree allele
  vector<int> ped_alleles(which_cc.size(), -1);
  for (unsigned int i = 0; i < ped_alleles.size(); i++)
    {
      int cc = which_cc[i];
      ped_alleles[i] = set_alleles[2*i+sampled_cc[cc]];
      if (ped_alleles[i] == 0)
	ped_alleles[i] = 1;
    }

  return ped_alleles;
}





// pick CC alleles uniformly from the 1-2 options
// if any allele possible, give it the major allele (i.e. 1)
//
// sample anything, but compute the rejection based probability
vector<int> IBDSetForest::SampleInfiniteAlleles()
{
  log_sample_prob = 0.0;



  // enumerate allele samples,
  // keep only the ones that are feasible (i.e. only a single 2 among all founders)

  vector<vector<int> > samples;



  // for each CC
  vector<int> sampled_cc(num_cc, -1);
  samples.push_back(sampled_cc);
  
  // assign cc's
  for (unsigned int i = 0;  i < which_cc.size();  i++)
    {
      int cc = which_cc[i];
      if (samples[0][cc] != -1)
	continue;

      vector<int> options;
      vector<int> alleles = GetAlleles(i);
      if (alleles[0] == 0 || alleles[1] == 0)
	{
	  options.push_back(0);
	}
      else // otherwise
	{
	  int num_alleles = 0;
	  for (unsigned int a = 0;  a < 2;  a++)
	    if (alleles[a] != -1)
	      num_alleles++;

	  // Check whether there are actually two different options
	  unsigned int diff_count = 0;
	  if (num_alleles == 2)
	    {
	      for (unsigned int j = 0;  j < which_cc.size();  j++)
		{
		  if (which_cc[j] == cc)
		    if (set_alleles[2*j+0] != set_alleles[2*j+1])
		      {
			diff_count++;
			break;
		      }
		}
	    }
	  unsigned int upper_bound = num_alleles;
	  if (num_alleles == 2)
	    upper_bound = 1 + diff_count;
	  for (unsigned int a = 0;  a < upper_bound;  a++)
	    options.push_back(a);
	}

      unsigned int num_curr_samples = samples.size();
      for (unsigned int j = 0;  j < num_curr_samples;  j++)
	{
	  //if (options.size() == 0)
	  //  samples[j][cc] = -1;
	  //else 
	  if (options.size() == 1)
	    samples[j][cc] = options[0];
	  else if (options.size() == 2)
	    {
	      samples[j][cc] = options[0];
	      vector<int> new_sample = samples[j];
	      new_sample[cc] = options[1];
	      samples.push_back(new_sample);
	    }
	}
    }



  // check each sample for feasibility
  // recording feasible sample indexes in a list
  vector<unsigned int> feasible_index;

  for (unsigned int j = 0;  j < samples.size();  j++)
    {
      //printf(" SAMPLE %i] ", j);
      // for each founder, count num of 2-alleles
      int two_alleles = 0;
      for (unsigned int i = 0;  i < which_cc.size(); i++)
	{
	  int cc = which_cc[i];

	  if (set_alleles[2*i+samples[j][cc]] == -1)
	  {
	    two_alleles = 99;
	    break;
	  }
	  
	  // check founder status
	  if (IsRoot(i)){
	    //printf("a%i,v%i, ", i, set_alleles[2*i+samples[j][cc]]);
	    if (2 == set_alleles[2*i+samples[j][cc]])
	      two_alleles++;
	  }
	}
      if (two_alleles <= 1)
	{
	  feasible_index.push_back(j);
	  //printf("   FEASIBLE");
	}
      //printf("\n");
    }


  // set probability for this sample
  log_sample_prob = -log(feasible_index.size());

  // choose randomly from the feasible indices
  int idx = -1;
  if (feasible_index.size() > 0)
    idx = (int) floor(drand48() * 324003043) % feasible_index.size();
  
  //printf("  choose %i with prob: 1/%i\n", idx, feasible_index.size());



  // copy selected alleles into an allele vector, one entry per pedigree allele
  vector<int> ped_alleles(which_cc.size(), -1);
  for (unsigned int i = 0; i < ped_alleles.size(); i++)
    {
      int cc = which_cc[i];
      if (idx == -1)
	ped_alleles[i] = -1;
      else
	{
	  int s = feasible_index[idx];
	  ped_alleles[i] = set_alleles[2*i+samples[s][cc]];
	  if (ped_alleles[i] == 0)
	    ped_alleles[i] = 1;
	}
    }


  return ped_alleles;
}








#endif
