

#ifndef FAMILYALG
#define FAMILYALG

#include <string.h>
#include "pedigree.h"



class Family;



// Had trouble with making an abstract class (i.e. it wouldn't compile)
class CompositeAlgorithm
{
 public:
  virtual void run(Family& family){};
  virtual void childRun(Family& family){};
  virtual void printAlgorithmNames(){printf("CA::pAN()\n");};
  virtual void setRange(unsigned int _start_marker, unsigned int _upper_bound_marker){};

  //CompositeAlgorithm(){printf("WARNING: Instantiating CompositeAlgorithm() could lead to unexpected behavior.\n");};
  //void run(Family& family){printf("WARNING: Running CompositeAlgorithm.run() could lead to unexpected behavior.\n"};
  //void childRun(Family& family){printf("WARNING: Running CompositeAlgorithm.childRun() could lead to unexpected behavior.\n");};
  //void printAlgorithmNames(){printf("WARNING: Running CompositeAlgorithm.printAlgorithNames() could lead to unexpected behavior.\n");};
  //void setRange(unsigned int _start_marker, unsigned int _upper_bound_marker){printf("WARNING: Running CompositeAlgorithm.setRange() could lead to unexpected behavior.\n");};
};




class FamilyAlgorithm : public CompositeAlgorithm
{
 public:
  FamilyAlgorithm(){
    backward_recursion = 0;
    child=NULL;
    strcpy(name, "Family Algorithm (generic)");
  };


  virtual void run(Family& family); 
  virtual void childRun(Family& family); // needs to run children first
  virtual void addChildAlgorithm(FamilyAlgorithm& alg);


  virtual void setRange(unsigned int _start_marker, unsigned int _upper_bound_marker);
  virtual void setForward();
  virtual void setBackward();
  virtual void printAlgorithmNames();


  bool isVerbose(){return verbose;}
  void setVerbose(bool v){verbose = v;}


 protected:
  char name[256]; // algorithm name
  bool verbose;

  unsigned int start_marker;
  unsigned int upper_bound_marker;
  int backward_recursion;

  CompositeAlgorithm* child;

};





void FamilyAlgorithm::setRange(unsigned int _start_marker, unsigned int _upper_bound_marker)
{
  //printf("[FamAlg] %s.setRange()\n", name);

  if (_start_marker <= _upper_bound_marker){
    start_marker = _start_marker;
    upper_bound_marker = _upper_bound_marker;
  } else {
    printf("FATAL ERROR in FamilyAlgorithm::setRange(%u,%u): markers out of bounds\n", _start_marker, _upper_bound_marker);
    exit(-1);
  }
}


void FamilyAlgorithm::setForward(){backward_recursion = 0;};
void FamilyAlgorithm::setBackward(){backward_recursion = 1;};



// runs children first
void FamilyAlgorithm::run(Family& family)
{
  //printf("[FamAlg] %s.run()\n", name);

  this->childRun(family);

  // block clean-up details...
}

void FamilyAlgorithm::childRun(Family& family)
{
  //printf("[FamAlg] %s.childRun()\n", name);

  if (child != NULL){
    child->setRange(start_marker, upper_bound_marker);
    child->childRun(family);
  }

  // alg details....
}



void FamilyAlgorithm::addChildAlgorithm(FamilyAlgorithm& alg)
{
  child = &alg;
}


void FamilyAlgorithm::printAlgorithmNames()
{
  if (child != NULL)
    {
      child->printAlgorithmNames();
    }
  printf("%s\n", name);
}


#endif
