
#ifndef STATES
#define STATES

#include <vector>
#include <algorithm>
#include <string>
#include <set>
#include <map>


using namespace std;


//
// List Vector
//   inheritance path indices
//   log probability entries
//
// Converts inheritance path vectors into integers
// Hamming distance btw indices is computed:
//    x = i (XOR) j
//    h = sum_k  x (XOR) e_k
// where e_k is a basis vector/integer (i.e. 10000, or 010000)
//

typedef map<long,double>::iterator GenoProbIterator;

class GenoProb
{
 public:
  GenoProb(int _num_nonfounders){
    num_nonfounders = _num_nonfounders;
    num_basis = (int) 2*num_nonfounders;
  };
  ~GenoProb(){
    vector.clear();
  };

  GenoProb(const GenoProb& g){
    num_nonfounders = g.num_nonfounders;
    num_basis = g.num_basis;
    vector = g.vector;
  };
  const GenoProb& operator= (const GenoProb& g){
    GenoProb* n = new GenoProb(g.num_nonfounders);
    n->vector = g.vector;
    return *n;
  };



  int size(){ return vector.size(); };

  GenoProbIterator begin(){ return vector.begin(); };
  GenoProbIterator end(){ return vector.end(); };

  long get_index(GenoProbIterator it){ return it->first; };
  double get_prob(GenoProbIterator it){ return it->second; };
  double get_prob(long index){ return vector[index]; };

  void set_prob(int index, double prob){ vector[index] = prob; };

  int hamming_dist(int index1, int index2);

  int get_basis(){ return num_basis; };


 protected:
  map<long, double> vector;

  int num_nonfounders;
  int num_basis;

};


int GenoProb::hamming_dist(int index1, int index2)
{
  int xor_value = index1 ^ index2;

  //printf("    xor %i + %i = %i\n", index1, index2, xor_value);

  int dist = 0;
  int base = 1;
  //printf("basis sum: ");
  for (int i = 0;  i < num_basis;  i++)
    {
      int mask = xor_value & base;
      if (mask > 0){
	//printf(" + 1 (base %i, mask %i)  ", base, mask);
	dist += 1;
      }
      base *= 2;
    }
  //printf("\n");

  //printf("hamming %i to %i = %i\n", index1, index2, dist);

  return dist;
}



// 
// List Vector
//   haplotype parental origin indices
//   entries are genomatrices
//   record of haplotypes' parental origins
//
// Haplo parental origin is determined by the allele sets for each 
// CC of the inheritance graph.  If the geno is heterozygous at that locus
// and if the inheritance graph dictates the orientation of the alleles, then
// the parental origin is fixed, otherwise, it is not fixed.
// We are only interested in scoring (summing probabilties) for indices with 
// the same haplotype origins.  
//
// probability indices: parental origin x inheritance paths
//

typedef map<long,GenoProb*>::iterator HaploProbIterator;

class HaploProb
{
 public:
  HaploProb(int _num_nonfounders){
    num_nonfounders = _num_nonfounders;
    num_basis = 2 * num_nonfounders;
  }

  HaploProb(const HaploProb& h){
    num_nonfounders = h.num_nonfounders;
    num_basis = h.num_basis;
    orientation_vector = h.orientation_vector;
  };
  const HaploProb& operator= (const HaploProb& h){
    HaploProb* n = new HaploProb(h.num_nonfounders);
    n->orientation_vector = h.orientation_vector;
    return *n;
  };


  int size(){ return orientation_vector.size(); };
  HaploProbIterator begin() { return orientation_vector.begin(); };
  HaploProbIterator end() { return orientation_vector.end(); };

  long get_orientation(HaploProbIterator it){ return it->first; };
  GenoProb* get_genoprob(HaploProbIterator it){ return it->second; };
  GenoProb* get_genoprob(int orientation){ 
    HaploProbIterator h = orientation_vector.find(orientation);
    if (h != orientation_vector.end())
      return h->second;
    else
      {
	//printf("RETURNED BLANK\n");
	GenoProb* g = new GenoProb(num_nonfounders);
	//orientation_vector[orientation] = g;
	return g;
      }
  };

  void set_genoprob(int orientation_index, GenoProb* g){ orientation_vector[orientation_index] = g; };


  double get_basis(){ return num_basis; };



 protected:
  map<long, GenoProb*> orientation_vector;

  int num_nonfounders;
  int num_basis;

};



#endif
