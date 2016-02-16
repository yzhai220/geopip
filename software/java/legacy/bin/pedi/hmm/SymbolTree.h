
/////////////////////////////////////
//
// Author: Javier Rosa
//   Date: Aug. 2007
//
//  Desc.: Maintains a symbol tree for boolean 
//         expressions.  Converts the symbol to 
//         tree to DNF on request, and provides 
//         an iterator to access all the clauses 
//         in the DNF expression.
//
//         Also provides some functionality
//         specific to haplotype boolean 
//         expressions, specifically:
//              UNKNOWN_LITERAL
//              cull()
//
//         To properlly use shuffle(), initialize 
//         the random number generator in the 
//         calling application:
//              srand48(time(NULL));
//
//
/////////////////////////////////////



#ifndef SYMBOLTREE_H
#define SYMBOLTREE_H
#define StartD 
#include <list>
#include <vector>

#define JRDEBUG 1
//REMOVE THIS LINE TO ENABLE DEBUG
#undef JRDEBUG

#ifdef JRDEBUG
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


#define UNKNOWN_LITERAL -33 // this litteral can take any value


using namespace std;
//Forward Declarations
class SymNode;
class SymTree;

// Enumerated type for Node Type 
// BTerm is the term that the node represents ie terms in (1^2)\/78 are 1,2,78
enum NType {BTERM = 0,SAND,SOR,SNOT,FTYPE,TTYPE,ETYPE};
// FTYPE == false and TTYPE == true
const int MODIFIED = -1;
const int READY = 1;
const int NOMORE = 0;


#define DUMMY -2
typedef long  int BTerm; // literal
const int NOPARENT = -1;
#define ChildList list<SymNode *>
#define CLIter ChildList::iterator

//Prototypes
//These are just a quick way to create a tree from binary ands
inline SymTree *and2(BTerm, BTerm);
inline SymTree *or2(BTerm, BTerm);
inline SymTree *and2(SymTree *, SymTree *);
inline SymTree *or2(SymTree *, SymTree *);
inline SymTree *or2(SymTree *a, BTerm j);
inline SymTree *and2(SymTree *a, BTerm j);
class SymNode 
{
  friend class SymTree;
 private: 
  void removeParens();
 public:
  //Type determines if it is an operator or a term
  NType t;
  BTerm bt;
  int index; //This is the index of this node in it's parent -1 otherwise
  ChildList children;

  //Constructor
  explicit SymNode(NType, BTerm);
  explicit SymNode(SymNode*);

  void addChild(SymNode*);
  void addChildren(SymNode*);

  //NOT FINISHED
  void removeChild(int);

  void destroyTree();
  void destroyNodeShallow(){
    this->children.clear();
    delete this;
  }

  SymNode *copyTree();
  bool cullNode();
  bool cullNode2();
  SymNode *toDNF();
  void printTree(int);
  int distribute();
  bool equals(SymNode *);
  friend SymNode *distr(SymNode *, SymNode *);
  
};
                                            
class SymTree			
{	
 private:
  SymTree *metaOp(SymTree *, NType);
  SymNode *root;		
  CLIter curAnd;
  int ready;
  void destroyTree();

 public:
  //Constructors			
  explicit SymTree(BTerm);
  explicit SymTree();
  //Copy constructor
  explicit SymTree(SymTree *);
  //Destructortron
  ~SymTree();

  //The argument trees are cleared and replaced with empty Trees
  SymTree *sand( SymTree *);	
  SymTree *sor( SymTree *);	 
  SymTree *sand( BTerm );	
  SymTree *sor( BTerm );
  void toDNF();	

  void cull();
  void enforceUnknownPrecedence();
  bool isEmpty();
  void printTree();
  friend void printNode(SymNode *,int);
  friend SymTree *shallowCopyAndClear(SymTree *);

  int get_next_pair(int&, int&);
  void rewind();
  void shuffle(); // randomly sort the children of the tree
};                              


//These are UNSTABLE use the returned values


inline SymTree *bkAnd2(SymTree *a, BTerm j){
  SymTree * result = a->sand(j);
  delete a;
  return result;
}

inline SymTree *bkOr2(SymTree *a, BTerm j){
  SymTree * result = a->sor(j);
  delete a;
  return result;
}
inline SymTree* bkAnd2(SymTree *a, SymTree *b){
  SymTree * result = a->sand(b);
  delete a;
  //delete b;
  return result;
}

inline SymTree* bkOr2(SymTree *a, SymTree *b){
  SymTree * result = a->sor(b);
  delete a;
  //delete b;
  return result;
}


//These are UNSTABLE use the returned values
inline SymTree *and2(SymTree *a, BTerm j){
  SymTree * result = a->sand(j);
  delete a;
  return result;
}

inline SymTree *or2(SymTree *a, BTerm j){
  SymTree * result = a->sor(j);
  delete a;
  return result;
}
inline SymTree* and2(SymTree *a, SymTree *b){
  SymTree * result = a->sand(b);
  delete a;
  delete b;
  return result;
}

inline SymTree* or2(SymTree *a, SymTree *b){
  SymTree * result = a->sor(b);
  delete a;
  delete b;
  return result;
}

//STABLE
inline SymTree * and2(BTerm i, BTerm j){
  //Should Validate BTerms
  SymTree *a = new SymTree(i);
  SymTree *b = new SymTree(j);
  SymTree *c = and2(a,b);
  return c;
}

inline SymTree * or2(BTerm i, BTerm j){
  //Should Validate BTerms
  SymTree *a = new SymTree(i);
  SymTree *b = new SymTree(j);
  SymTree *c = or2(a,b);
  return c;
}


#endif
