
#ifndef IBDALG
#define IBDALG

#include "pedigree.h"
#include "familyalgorithm.h"
#include "ibd_set_forest.h"
#include "state_space.h"


#define IBDALGDEBUG 1
//REMOVE THIS LINE TO ENABLE DEBUG
#undef IBDALGDEBUG

#ifdef IBDALGDEBUG
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





// Prevents underflow by being large enough,
// but might induce overflow if it is too large
#define MULTIPLIER log(64);
double EXPMULT = 64.0;




class IBDAlg : public FamilyAlgorithm
{
 public:
  IBDAlg()
  {
    backward_recursion = 0;
    child = NULL;
    rec_rate = 0;
    strcpy(name, "IBD Graph Algorithm");
  };
  virtual ~IBDAlg()
  {
  };


  virtual void run(Family& family);
  virtual void childRun(Family& family);

  void setForward(){backward_recursion = 0;};
  void setBackward(){backward_recursion = 1;};



 protected:

  void enumerate_paths(Family& family,
		       vector<int>& partial_path, 
		       const unsigned int order_index);
  void baseCase(Family& family, 
		vector<int>& partial_path);

  void marginalTransitionProb(Family& family);
  void marginalProb(Family& family, int two_nf);


  // one entry per untyped person, with a map containing sums of likelihoods
  vector<map<string, double> > marginals;  
  vector<map<string, double> > marginal_entropy;  

  IBDSetForest* forest;

  int two_nf;
  int nf;
  double rec_rate;
  unsigned int num_snps;

};












void IBDAlg::run(Family& family)
{
  this->childRun(family);

  // clean up any loci variables, when this alg is run
  // as its own parent
}




void IBDAlg::childRun(Family& family)
{
  //printf("IBDAlg::childRun()\n");

  // first run child
  if (child != NULL)
    {
      child->setRange(start_marker, upper_bound_marker);
      child->childRun(family);
    }


  // Enumerate the inheritance paths in reverse topological order.
  // Keep an updated IBD graph for each partial enumeration
  // (i.e. an geno-labeled edge between inheritance paths 
  // whenever they intersect in a genotyped person).
  // Bound the enumeration of inheritance paths if the IBD graph 
  // contains no legal allele states.


  //#ifdef DEBUG
  //printf("  ibd alg (fam %i, start_marker %i, backward_recursion %i)\n", family.family_id, start_marker, backward_recursion);
  //fflush(stdout);
  //#endif


  // For haplotype data give the hap ordering determined by 
  // the heterozygous loci.


  if (start_marker+1 !=  upper_bound_marker)
    {
      printf("FATAL ERROR: IBDAlg given range for multi-marker block... should be single marker\n");
      exit(-36787);
    }



  // class variables specific to this family
  if (start_marker != 0)
    {
      rec_rate = family.pedigree->recombination_rates[start_marker-1];
    }
  if (verbose)
    printf("\n\nFam %i, REC RATE: %f\n\n", family.family_id, rec_rate);

  nf = family.nonfounder.size();
  two_nf = 2*family.nonfounder.size();
  num_snps = family.members[0]->genotypes.size();


  family.curr_inheritance.clear();
  family.curr_inheritance.resize(pow(2,two_nf));

  //family.inheritance.clear();
  //family.inheritance.resize(num_snps, vector<IBDSetForest>(pow(2,two_nf)));
  family.marginal_indicators.clear();
  family.marginal_indicators.resize(num_snps, vector<double>());


  if (backward_recursion)
    {
      if (start_marker+1 == num_snps)
	{
	  family.geno_entropy = 0;
	}

      if (family.emission_prob != NULL)
	{
	  delete family.emission_prob;
	  family.emission_prob = NULL;
	}
      family.emission_prob = new GenoProb(nf);

      if (family.geno_backward_current != NULL)
	{
	  delete family.geno_backward_current;
	  family.geno_backward_current = NULL;
	}
      if (start_marker+1 == num_snps && family.geno_backward != NULL)
	{
	  delete family.geno_backward;
	  family.geno_backward = NULL;
	}


      family.geno_backward_current = family.geno_backward; // current marker
      if (start_marker != 0)
      {
	  family.geno_backward = new GenoProb(nf); // previous marker, i.e. in the genome sequence
      }
    }
  else
    {
      GenoProb* this_geno_marker = new GenoProb(nf);
      family.geno_dynprog.push_back(this_geno_marker);
    }



  // begin recursive enumeration
  forest = new IBDSetForest(family.members.size(), 2*family.founder.size());
  for(unsigned int i = 0;  i < family.members.size(); i++){
    int j = family.members[i]->order+1;
    int allele1 = family.members[i]->genotypes[start_marker].allele1;
    int allele2 = family.members[i]->genotypes[start_marker].allele2;
    int num_children = family.members[i]->children.size();
    forest->InsertNode(j, i, allele1, allele2, num_children);
  }
  forest->SetNumAlleles(2);
  vector<int> partial_path(two_nf, -1);
  enumerate_paths(family, partial_path, 0);


  // REPORTING FOR RESULTS OF Recursive Step
  if (!backward_recursion)
    {
      // FORWARD probability
      //if (start_marker+1 == num_snps)
	{
	  // sum the last marker's probabilities to obtain 
	  // the prob of the observed data
	  //GenoProb* last_marker = family.geno_dynprog[family.geno_dynprog.size()-1];
	  GenoProb* last_marker = family.geno_dynprog[start_marker];
	  double prob = 0;
	  double sum = 0;
	  //printf("last_marker.size() = %i\n", last_marker->size());
	  for (GenoProbIterator prev_iterator = last_marker->begin();  
	       prev_iterator != last_marker->end();  
	       prev_iterator++)
	    {
	      double index_prob = last_marker->get_prob(prev_iterator);
	      double index_prob_normalized = index_prob/EXPMULT/last_marker->size();
	      prob += index_prob / last_marker->size();
	      //int index = last_marker->get_index(prev_iterator);
	      //printf("   %i has prob %f\n", index, index_prob_normalized);
	      sum += index_prob_normalized;
	    }
	  //printf("      sum: %f\n", sum);
	  family.geno_log_prob_data = log(prob) - num_snps * MULTIPLIER;
	  

	  if (verbose)
	    printf("\n\nfamily %i] Probability of  geno observation: %f\n", family.family_id, family.geno_log_prob_data);



	}
    }
  else
    {
      // BACKWARD probability
      
      // compute the local transition probability
      if (start_marker != 0)
	marginalTransitionProb(family);

      // compute the marginals
      marginalProb(family, two_nf);
    }


  // think there is memory leak unless we delete
  delete forest;

  // fixes a memory error
  if (start_marker == 0 && family.geno_backward != NULL)
    {
      delete family.geno_backward;
      family.geno_backward = NULL;
    }

  if (start_marker == 0 && family.emission_prob != NULL)
    {
      delete family.emission_prob;
      family.emission_prob = NULL;
    }
  
  //printf("\n\n");
}




void IBDAlg::marginalTransitionProb(Family& family)
{
  ////////////
  // Geno


  vector<double> recomb_prob(two_nf+1,0);


  // deal with end effect
  GenoProb* bw_current = family.geno_backward_current;
  int isVectorOnes = false;
  if (start_marker+1 == num_snps)
    {
      isVectorOnes = true;
      // because we will need to iterate over the feasible inheritance paths
      bw_current = family.geno_dynprog[start_marker];
    }



  double sum = 0;
  GenoProb* prev_marker = family.geno_dynprog[start_marker-1];  
  for (GenoProbIterator prev_iterator = prev_marker->begin();  
       prev_iterator != prev_marker->end();  
       prev_iterator++)
    {
      for (GenoProbIterator curr_iterator = bw_current->begin();
	   curr_iterator != bw_current->end();
	   curr_iterator++)
	{
	      double curr_prob = bw_current->get_prob(curr_iterator);
	      // FIX THIS BAD IMPLEMENTATION (lack of locality)
	      // should have values in bw_current even for last SNP
	      if (isVectorOnes)
		curr_prob = 1*EXPMULT;
	      int curr_path = bw_current->get_index(curr_iterator);
	      double emission_prob = family.emission_prob->get_prob(curr_path);

	      double prev_prob = prev_marker->get_prob(prev_iterator);
	      int prev_path = prev_marker->get_index(prev_iterator);

	      int h = prev_marker->hamming_dist(prev_path, curr_path);


	      double log_prob = log(pow(rec_rate,h)) + log(pow(1-rec_rate,two_nf-h));
	      log_prob += log(prev_prob) + log(curr_prob) + emission_prob;
	      log_prob -= num_snps * MULTIPLIER; // all the correction that is required
	      log_prob -= family.geno_log_prob_data;

	      double prob = exp(log_prob);
	      sum += prob;

	      recomb_prob[h] += prob;

	      //printf("    trans (%i,%i) has E = %f\n", prev_path, curr_path, prob);
	}
    }
  //printf("      sum %f\n", sum);


  double entropy = 0;
  
  family.geno_entropy += entropy;
}


void IBDAlg::marginalProb(Family& family, int two_nf)
{
  ////////////
  // Geno


  // deal with end effect
  GenoProb* bw_current = family.geno_backward_current;
  int isVectorOnes = false;
  if (start_marker+1 == num_snps)
    {
      isVectorOnes = true;
      // because we will need to iterate over the feasible inheritance paths
      bw_current = family.geno_dynprog[start_marker];
    }


  double sum = 0;
  GenoProb* prev_marker = family.geno_dynprog[start_marker];

  // vector for marginal
  vector<double> marginal(pow(2, prev_marker->get_basis()),0);

  vector<double> marginal_indicator(4*two_nf,0);


  GenoProbIterator prev_iterator;
  GenoProbIterator curr_iterator;
  for (prev_iterator = prev_marker->begin(),
	 curr_iterator = bw_current->begin();
       prev_iterator != prev_marker->end(),
	 curr_iterator != bw_current->end();
       prev_iterator++,
	 curr_iterator++)
    {
      double curr_prob = bw_current->get_prob(curr_iterator);
      // FIX THIS BAD IMPLEMENTATION (lack of locality)
      // should have values in bw_current even for last SNP
      if (isVectorOnes)
	curr_prob = 1*EXPMULT;
      int curr_path = bw_current->get_index(curr_iterator);
      //double emission_prob = family.emission_prob->get_prob(curr_path);
      
      double prev_prob = prev_marker->get_prob(prev_iterator);
      //int prev_path = prev_marker->get_index(prev_iterator);
      
      double log_prob = 0;
      log_prob += log(prev_prob) + log(curr_prob); 
      log_prob -= (num_snps +1)* MULTIPLIER; // all the correction that is required
      log_prob -= family.geno_log_prob_data;
      
      double prob = exp(log_prob);
      sum += prob;

      // full inheritance path posterior
      if ((unsigned) curr_path < marginal.size())
	marginal[curr_path] += prob;

      //printf("    state (%i) has prob = %f\n", curr_path, prob);      




      // inheritance indicator posterior marginal
      // loop over flat_set_forest, 
      // for each inheritance pointer, sum prob with it indicator index in vector
      IBDSetForest *p = &family.curr_inheritance[curr_path];
      for (int i = 0; i < two_nf; i++)
	{
	  for (int parent = 0 ;  parent  < 2;  parent++)
	    {
	      marginal_indicator[4*i +2*parent +p->GetAlleleIndicatorPed(i,parent)] += prob;
	    }
	}



    }
  //printf("      sum %f\n", sum);


  // normalize all the distributions

  for (unsigned int h = 0;  h < marginal.size();  h++)
    {
      marginal[h] /= sum;
      //printf("    state (%i) has prob = %f\n", h, marginal[h]);
    }

  for (unsigned int i = 0 ;  i < marginal_indicator.size();  i++)
    {
      marginal_indicator[i] /= sum;
    }


  // store data in vectors
  //family.inheritance[start_marker] = family.curr_inheritance;
  family.marginal_indicators[start_marker] = marginal_indicator;


  double entropy = 0;
  
  family.geno_entropy += entropy;
}


void IBDAlg::baseCase(Family& family, 
		      vector<int>& partial_path)
{

  // emission probability for observed data given this inheritance path
  // and the input founder allele frequencies
  double emission_prob = forest->GetLogEmissionProb(); 


  double emission_weight = 0;
  if (start_marker < family.indicator_weights.size())
    emission_weight = forest->GetLogEmissionWeight(family.indicator_weights[start_marker]);

  emission_prob += emission_weight;
  


  if (emission_prob <= -999999999)
    {
      printf("ERROR: baseCase should not have been reached, but valid allele set: %i.\n",  forest->validAlleleSets());
      exit(-5435634);
    }



  // Convert the partial path into an integer.
  int path_index = 0;
  int base = 1;
  //printf("IBDAlg::baseCase()   likelihood %f\n", exp(log_likelihood));
  for (unsigned int i = 0;  i < partial_path.size(); i++)
  {
    deboutnn(partial_path[i]<<" ");
    
    path_index += partial_path[i] * base;
    base *= 2;
  }
  deboutln(" => "<<path_index<<" with prob "<<exp(emission_prob));



  // save the current inheritance indicators
  family.curr_inheritance[path_index] = *forest;




  ///////////
  //  GENO
  ///////////
  // Loop over the indices of the previous marker
  // and do sum-product to compute the probability for 
  // this index's probability.


  if (!backward_recursion)
    {
      //
      // FORWARD ITERATION
      //

      if (family.geno_dynprog.size() <= start_marker)
	{
	  printf("ERROR IN creating dynprog array: size %i, start_marker %i\n", (int) family.geno_dynprog.size(), start_marker);
	  exit(-43534);
	}


      GenoProb* this_marker = family.geno_dynprog[start_marker];
      if (start_marker == 0)
	{
	  // separately set the first marker with equal prob prior
	  //double prob = log_likelihood + emission_prob;
	  double prob = emission_prob + MULTIPLIER;
	  this_marker->set_prob(path_index, exp(prob));
	}
      else
	{
	  double prob = 0;
	  GenoProb* prev_marker = family.geno_dynprog[start_marker-1];
	  for (GenoProbIterator prev_iterator = prev_marker->begin();  
	       prev_iterator != prev_marker->end();  
	       prev_iterator++)
	    {
	      int prev_path = prev_marker->get_index(prev_iterator);
	      double prev_log_prob = log(prev_marker->get_prob(prev_iterator));
	      int h = prev_marker->hamming_dist(prev_path, path_index);

	      //printf("       geno path %i, pprob %f, dist %i, r %f\n", prev_path, prev_log_prob, h, rec_rate);
	      double term = log(pow(rec_rate,h)) + log(pow(1-rec_rate,two_nf-h)) + prev_log_prob + emission_prob + MULTIPLIER;
	      //term -= log_likelihood;
	      prob += exp(term);
	    }
	  //deboutln("  geno path "<<path_index<<", prob "<< prob);
	  this_marker->set_prob(path_index, prob);
	}
    }
  else
    {
      // 
      // BACKWARD ITERATION
      //

      //printf("start_marker %i    dp size: %i\n", start_marker, family.geno_dynprog.size());

      if (start_marker != 0)
	{
	  // marker order:  geno_backward, geno_backward_current
	  // we sum over elements of geno_backward_current to obtain value for geno_backward

	  // emission prob for geno_backward_current which has prob. for start_marker's locus
	  family.emission_prob->set_prob(path_index, emission_prob);

	  double backward_log_prob = 0 + MULTIPLIER; // vector of 1's for final marker
	  if (start_marker+1 != num_snps)
	    backward_log_prob = log(family.geno_backward_current->get_prob(path_index));
	  //else
	  //family.geno_backward_current->set_prob(path_index, 1);

	  GenoProb* prev_marker = family.geno_dynprog[start_marker-1];
	  for (GenoProbIterator prev_iterator = prev_marker->begin();
	       prev_iterator != prev_marker->end();  
	       prev_iterator++)
	    {
	      int prev_path = prev_marker->get_index(prev_iterator);
	      int h = prev_marker->hamming_dist(prev_path, path_index);

	      // sum this term into the geno_backward prob for this prev_path
	      double term = 0;
	      term += emission_prob + backward_log_prob;
	      term += log(pow(rec_rate,h)) + log(pow(1-rec_rate,two_nf-h));
	      term += MULTIPLIER;

	      double prev_prob = family.geno_backward->get_prob(prev_path);
	      double prob = prev_prob + exp(term);
	      family.geno_backward->set_prob(prev_path,prob);

	      //printf("       path %i, pprob %f, dist %i, r %f\n", prev_path, prev_prob, h, rec_rate);
	      //printf("    prob for %i is %f\n", prev_path, prob);
	    }
	} // end if (start_marker != 0)

    }



   
}




void IBDAlg::enumerate_paths(Family& family,
			     vector<int>& partial_path, 
			     const unsigned int order_index)
{
  
  //#ifdef DEBUG
  //printf("  enumerate_paths (fam %i, order %i, log_likelihood %f), start_marker %i, backward_recursion %i\n", family.family_id, order_index, log_likelihood, start_marker, backward_recursion);
  //fflush(stdout);
  //#endif






  // base case: we have reached the top of the pedigree
  if (order_index >= family.order.size())
    {
      //printf("  IBDAlg::enumerate_paths() BASE w/ ll = %f\n", log_likelihood);
      bool isCompatible = forest->validAlleleSets();
      if(isCompatible){
	baseCase(family, partial_path);
      }

      return;
    }

 

  // for this order_index (in bottom_up order)
  OrderItem order_item = family.order[order_index];
  Individual* person = order_item.individual;
  Individual* mother = person->mother;
  Individual* father = person->father;


  //#ifdef DEBUG
  //    printf("    (ibd) Processing indiv %i:%i (", family.family_id, order_item.individual_id);
  //    if (father != NULL && mother != NULL)
  //	printf("NONfounder)\n");
  //    else
  //    printf("founder)\n");
  //    fflush(stdout);
  //#endif



  // case: founder
  if (mother == NULL || father == NULL)
    { 
      enumerate_paths(family, partial_path, order_index+1);
    }
  // case: non-founder
  else
    {
      int new_index = order_index+1;
      bool isCompatible = true;
      
      partial_path[2*order_index] = 0;
      partial_path[2*order_index+1] = 0;
      //printf("      person i%i o%i  father i%i,o%i:%i  mother i%i,o%i:%i\n", person->individual_id, person->order+1, father->individual_id, father->order+1, partial_path[2*order_index], mother->individual_id, mother->order+1, partial_path[2*order_index+1]);
      isCompatible = forest->Union(person->order+1, father->order+1, partial_path[2*order_index], mother->order+1, partial_path[2*order_index+1]);
      if(isCompatible){
	enumerate_paths(family, partial_path, new_index);
	}
      forest->ReverseUnion(person->order+1);

      partial_path[2*order_index] = 0;
      partial_path[2*order_index+1] = 1;
      //printf("      person i%i o%i  father i%i,o%i:%i  mother i%i,o%i:%i\n", person->individual_id, person->order+1, father->individual_id, father->order+1, partial_path[2*order_index], mother->individual_id, mother->order+1, partial_path[2*order_index+1]);
      isCompatible = forest->Union(person->order+1, father->order+1, partial_path[2*order_index], mother->order+1, partial_path[2*order_index+1]);
      if(isCompatible){
	enumerate_paths(family, partial_path, new_index);
	}
      forest->ReverseUnion(person->order+1);


      partial_path[2*order_index] = 1;
      partial_path[2*order_index+1] = 0;
      //printf("      person i%i o%i  father i%i,o%i:%i  mother i%i,o%i:%i\n", person->individual_id, person->order+1, father->individual_id, father->order+1, partial_path[2*order_index], mother->individual_id, mother->order+1, partial_path[2*order_index+1]);
      isCompatible = forest->Union(person->order+1, father->order+1, partial_path[2*order_index], mother->order+1, partial_path[2*order_index+1]);
      if(isCompatible){
	enumerate_paths(family, partial_path, new_index);
	}
      forest->ReverseUnion(person->order+1);

      
      partial_path[2*order_index] = 1;
      partial_path[2*order_index+1] = 1;
      //printf("      person i%i o%i  father i%i,o%i:%i  mother i%i,o%i:%i\n", person->individual_id, person->order+1, father->individual_id, father->order+1, partial_path[2*order_index], mother->individual_id, mother->order+1, partial_path[2*order_index+1]);
      isCompatible = forest->Union(person->order+1, father->order+1, partial_path[2*order_index], mother->order+1, partial_path[2*order_index+1]);
      if(isCompatible){
	enumerate_paths(family, partial_path, new_index);
	}
      forest->ReverseUnion(person->order+1);

      // unwind
      partial_path[2*order_index] = -1;
      partial_path[2*order_index+1] = -1;
    }

}









#endif
