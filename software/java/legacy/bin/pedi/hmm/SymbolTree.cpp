#include "SymbolTree.h"
#include <iostream>
#include <algorithm>
#include <stdlib.h>
#include <stdio.h>

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

void printNode(SymNode * node, int level);

inline SymTree *shallowCopyAndClear(SymTree *tree);
void noDuplicateNodes(ChildList &);
inline SymNode *newSand();
inline SymNode *newSor();
inline SymNode *newEmpty();
void testIterator(SymTree *tree);
SymTree::SymTree( BTerm bt ){
  this->root = new SymNode(BTERM,bt);
  ready = -1;
}
SymTree::SymTree(SymTree* tree){
  this->root = tree->root->copyTree();
  ready = -1;
}
SymTree::SymTree(){
  root = new SymNode(ETYPE,-1);
  ready = -1;
}
SymTree::~SymTree(){
  this->root->destroyTree();
  //delete root;
  root = NULL;
}
void SymTree::toDNF(){
  ready = -1;
  if( !this->isEmpty() || this->root->t != BTERM){
    SymNode* dnf_result = this->root->toDNF();
    assert (dnf_result == this->root);  // BK check for memory consistency
    //noDuplicateNodes(this->root->children);
    //this->root = dnf_result;
  }
}		
bool SymTree::isEmpty(){
  assert(root!=0);
  return root->t == ETYPE;
}
/*
Clears the tree and returns a shallow copy of it

*/

SymTree *shallowCopyAndClear(SymTree *tree){
  SymTree *returnMe = new SymTree(tree->root->bt);
  returnMe->ready = MODIFIED;
  tree->ready = MODIFIED;
  returnMe->root->t = tree->root->t;
  returnMe->root->addChildren(tree->root);
  tree->root->children.clear();
  tree->root->destroyTree();
  tree->root = newEmpty();
  return returnMe;
}




inline SymNode *newSand(){
  return new SymNode(SAND,-1);
}
inline SymNode *newSor(){
  return new SymNode(SOR,-1);
}
inline SymNode *newEmpty(){
  return new SymNode(ETYPE,-1);
}

///////////////////////
// Post-condition (memory management):
//   Return a result which has pointers to whichever 
//         nodes from this and other that are retained.
//   Nodes this and other are destroyed where appropriate
//         and replaced with newEmpty() nodes.
//   
SymTree *SymTree::metaOp(SymTree *other,NType op){
  assert(op == SOR || op == SAND);
  other->ready = MODIFIED;
  this->ready = MODIFIED;

  if(this->isEmpty()){
    return shallowCopyAndClear(other);
  }else if(other->isEmpty()){
    return shallowCopyAndClear(this);
  }

  if( op == SAND && root->t == FTYPE){
    other->root->destroyTree();
    other->root = newEmpty();
    return shallowCopyAndClear(this);
  }else if(op == SAND && other->root->t == FTYPE){
    this->root->destroyTree();
    this->root = newEmpty();
    return shallowCopyAndClear(other);
  }


  NType thisT = this->root->t;
  NType otherT = other->root->t;
  SymNode *newNodePtr;
  //For simplicity make sure that both trees have operators for roots
  if(otherT == BTERM){
    newNodePtr = new SymNode( op, DUMMY );
    newNodePtr->addChild(other->root);
    other->root = newNodePtr;            
    otherT = op;
  }
  if(thisT == BTERM){
    newNodePtr = new SymNode( op, DUMMY );
    newNodePtr->addChild(this->root);
    this->root = newNodePtr;            
    thisT = op;
  }

  if(thisT == op && otherT == op){
    this->root->addChildren(other->root);

    other->root->children.clear();
    delete other->root;
    other->root = newEmpty();

    return shallowCopyAndClear(this);
  }else if(thisT != op && otherT != op ){
    newNodePtr = new SymNode(op, -1);
    newNodePtr->addChild(this->root);
    newNodePtr->addChild(other->root);

    this->root = newNodePtr;
    other->root = new SymNode(ETYPE,-1);

    return shallowCopyAndClear(this);
  }else if( thisT != op && otherT == op ){
    other->root->addChild(this->root);

    this->root = NULL;
    this->root = newEmpty();

    return shallowCopyAndClear(other);
  } else if(thisT == op && otherT != op){
    this->root->addChild(other->root);

    other->root = NULL;
    other->root = newEmpty();
    
    return shallowCopyAndClear(this);
  }
  cerr <<"error in metaOp\n"<< endl;
  assert(0);
  return 0;
  
}

SymTree * SymTree:: sor(SymTree * other)
{
  return this->metaOp(other,SOR);
}
SymTree * SymTree:: sor(BTerm bt){
  SymTree* tree_bt = new SymTree(bt);   // BK fixed
  SymTree* rtn = this->sor(tree_bt);    // a memory
  delete tree_bt;                       // error (8/6/2007)
  return rtn;
}
SymTree *SymTree::sand(SymTree * other){
  return this->metaOp(other,SAND);
}
SymTree *SymTree::sand(BTerm bt){
  SymTree* tree_bt = new SymTree(bt);   // BK fixed
  SymTree* rtn = this->sand(tree_bt);   // a memory
  delete tree_bt;                       // error (8/6/2007)
  return rtn;
}

void SymTree::printTree(){
  cout << "\nSTART TREE\n" << endl;
  if(this->root != 0){
    this->root->printTree(0);
  }else {
    cout << "root is null.(ie false)" << endl;
  }
  cout << "\nEND TREE\n" << endl;
}
void SymNode::printTree(int level){
  CLIter iter, end;
  printNode(this,level);
  iter = children.begin();
  end = children.end();
  while(iter != end){
    (*iter)->printTree(level + 1);
    iter++;
  }
}





void printNode(SymNode * node, int level){
  for(int i = 0; i < level; i++){
    cout << " ";
  } 
  if( node->t == SAND ){
    cout << "and";
  }else if( node->t == SOR ){
    cout << "or.";
  }else if( node->t == SNOT ){
    cout << "not";
  }else if(node->t == ETYPE){
    cout << "E";
  }else if(node->t == FTYPE){
    cout << "F";
  }else if(node->t == TTYPE){
    cout <<"T";
  }else if( node->t == BTERM ){
    cout << node->bt;
  }else{
    cout << "ndef";
  }
  cout << "\n" << flush;
}
//SymNode operations
//Constructor
SymNode::SymNode(NType t, BTerm bt) : children()
{
  //Should perform checking to ensure that the NType value is correct
  this->t = t;
  this->bt = bt;
  this->index = -1;
  return;
}

SymNode::SymNode(SymNode * node) : children()
{
  this->t = node->t;
  this->bt = node->bt;
  this->index = node->index;  
  return;
}

void SymTree::destroyTree(){
  assert(root != 0);
  ready = MODIFIED;
  root->destroyTree();
  //delete root;
}

void SymNode::destroyTree(){  
  while(children.size() > 0){
    (*children.begin())->destroyTree();
    //delete (*children.begin());
    children.pop_front();
  }
  delete this;
}


////////////////////////////////////////////////////////////
//
// Similar to cull, but now when we find a pair with an CU, 
// where C is any literal and U is the unknown, we remove
// all other pairs involving C.
//
// Pre-Condition:
//   1) Cull has been called.
//
// Post-Condition:
//   1) Invalidates any use of get_next_pair() that is called 
//      in code that wraps calls to this function.
//   2) Removes all node-pairs that contain C and not U, 
//      when CU appears in the list of pairs.
//
void SymTree::enforceUnknownPrecedence()
{
  if (ready == MODIFIED)
    {
      printf("FATAL ERROR: SymTree::cull() must precede SymTree::enforceUnknownPrecedence().");
      assert(0);
    }
  

  // Iterate through the list of pairs
  // to make a list of literals appearing with U
  vector<int> mates_of_unknown;

  this->rewind();
  int a, b;
  while (get_next_pair(a, b))
    {
      if (a == UNKNOWN_LITERAL && b != UNKNOWN_LITERAL)
	  mates_of_unknown.push_back(b);
      else if (a != UNKNOWN_LITERAL && b == UNKNOWN_LITERAL)
	  mates_of_unknown.push_back(a);
      else if (a == UNKNOWN_LITERAL && b == UNKNOWN_LITERAL)
	{ // both a and b are unknown
	  this->root->destroyTree();
	  this->root = new SymNode(BTERM,UNKNOWN_LITERAL);
	  return;
	}
    }

  if(root->t == SAND || root->t == BTERM){
    // do nothing, because a single pair has no precedence
  }else if(root->t == SOR){
    // check the children
    CLIter cur = root->children.begin(); 
    while( cur != root->children.end() )
      {
	int deleteNode = false;
	SymNode* andNode = (*cur);

	if (andNode->t == BTERM)
	  {
	    vector<int>::iterator found = find(mates_of_unknown.begin(), 
					       mates_of_unknown.end(), 
					       (int) andNode->bt);
	    if (found != mates_of_unknown.end())
	      { // this is equivalent to CU where C is value of child1
		deleteNode = true;
	      }
	  }	
	else
	  { 
	    assert(andNode->t == SAND);

	    // There must are exactly two children of this node
	    assert(andNode->children.size() == 2);
	    
	    // Enforce precedence on the children
	    SymNode* child1 = *(andNode->children.begin());
	    SymNode* child2 = *(++andNode->children.begin());

	    if (child1->bt != UNKNOWN_LITERAL 
		&& child2->bt != UNKNOWN_LITERAL)
	      {	    
		// find child1->bt in our list
		vector<int>::iterator found = find(mates_of_unknown.begin(), 
						   mates_of_unknown.end(), 
						   (int) child1->bt);
		if (found != mates_of_unknown.end())
		  { // this is equivalent to CU where C is value of child1
		    deleteNode = true;
		  }
		
		// check child2 
		if (deleteNode == false)
		  {
		    found = find(mates_of_unknown.begin(), 
				 mates_of_unknown.end(), child2->bt);
		    if (found != mates_of_unknown.end())
		      {
			deleteNode = true;
		      }
		  }
	      }
	  } // end else if SAND node

	if (deleteNode)
	  {
	    (*cur)->destroyTree();
	    cur = root->children.erase(cur);
	    cur--;
	  }
	
	cur++;
      } // end for each child

  } // end if root is a SOR

}



void SymTree::cull(){
  //FINISH ME
  //NEEDS TO BE IN DNF BEFORE CALLING
  debout("c1");
  assert(root != 0);
  if(this->isEmpty() || this->root->t == FTYPE || this->root->t == TTYPE){
    ready = 0;
    return;
  }

  if(root->t == SAND){
    debout("c2");
    if(root->cullNode()){
      ready = 1;
      return;
    } else{
      root = new SymNode(FTYPE,-1);
      ready = 0;
      return;
    }
  }else if(root->t == SOR){
    debout("c3");
    CLIter cur = root->children.begin(); 
    while( cur != root->children.end() ){
      assert((*cur)->t == BTERM || (*cur)->t == SAND);
      if( (*cur)->cullNode() ){
	cur++;
      }else{
	cur = root->children.erase(cur);
      }
    }

    //All the nodes have been culled or kept and simplified
    noDuplicateNodes(root->children);

    if(root->children.size() > 1){
      //We have multiple children
      debout("c4");
      ready = 1;
      rewind();
      return; 
    }else if(root->children.size() == 1){
      //Only one "and" left, so make it the root
      debout("c5");
      deb(root->printTree(0));
      SymNode *soleChild = *(root->children.begin());
      root->children.pop_front();
      root->addChildren(soleChild);
      soleChild->children.clear();
      root->t = soleChild->t;
      root->bt = soleChild->bt;
      soleChild->destroyTree();
      ready = 1;
      return;
    }else{
      //ROOT HAS NO CHILDREN so they've all been falsed
      debout("c6");
      root->destroyTree();
      root = new SymNode(FTYPE,-1);
      ready = 0;
      return;
    }
  }else if(root->t == BTERM){
    ready = 1;
    return;
  }
}


/*
  Returns true if the node is not culled and modifies it accordingly
  If it is culled it deletes the node and returns false to indicate that 
  pointers to the destroyed node should be deleted.

  Post-condition: 'and' clades with BTERM children are SORTED
                  while duplicate literals in the 'and' are removed
*/

bool SymNode::cullNode(){
  debout("cn1");
  if(t == BTERM){
    return true;
  }
  //Keep in mind that the tree should already be in DNF
  assert(t == SAND && children.size() > 0);
  noDuplicateNodes(children);
  if(children.size() == 1){
    debout("cn2");
    assert((*children.begin())->t == BTERM);
    this->t = (*children.begin())->t;
    this->bt = (*children.begin())->bt;
    (*children.begin())->destroyTree();
    children.clear();
    return true;
  }else if (children.size() == 3) {
    debout("cn3");
    // Modification by BK: inserted this case (8/1/2007)
    // if one of the children is UNKNOWN_LITERAL, then remove that literal
    int removedUnknown = false;
    CLIter it = children.begin();
    while(it != children.end()){
      assert((*it)->t == BTERM);
      CLIter e = it;
      it++;
      if ((*e)->bt == UNKNOWN_LITERAL){
	removedUnknown = true;
	(*e)->destroyTree();
	children.erase(e);
      }
    }
    if (!removedUnknown){
      this->destroyTree();
      return false;
    }
  }else if ( children.size()>3 ){
    // Note: all the children will be distinct, due to noDuplicateNodes()
    debout("cn4");
    this->destroyTree();
    return false;
  }//Children size == 2
  debout("cn4");

  // Sort the 2-tuple
  CLIter one,two;
  one = two = children.begin(); 
  two++;

  if( (*one)->bt > (*two)->bt ){
    BTerm temp = (*two)->bt; 
    (*two)->bt = (*one)->bt;
    (*one)->bt = temp;
  }
  return true;
}




/*
  Returns true if the node is not culled and modifies it accordingly
  If it is culled it deletes the node and returns false to indicate that 
  pointers to the destroyed node should be deleted.

  Post-condition: AND clades with BTERM children are SORTED
                  while duplicate literals in the AND are removed
                  An OR node who might have child 1 will have child (1 /\ 1).
                  This satisfies the assumptions for culling in distribute().
*/

bool SymNode::cullNode2(){
  debout("cn1");
  if(t == BTERM){
    return true;
  }
  //Keep in mind that the tree should already be in DNF
  assert(t == SAND && children.size() > 0);
  noDuplicateNodes(children);
  if(children.size() == 1){
    debout("cn2");
    assert((*children.begin())->t == BTERM);
    this->addChild(new SymNode(BTERM, ((*children.begin())->bt)));
    //this->t = (*children.begin())->t;
    //this->bt = (*children.begin())->bt;
    //(*children.begin())->destroyTree();
    //children.clear();
    return true;
  }else if (children.size() == 3) {
    debout("cn3");
    // Modification by BK: inserted this case (8/1/2007)
    // if one of the children is UNKNOWN_LITERAL, then remove that literal
    int removedUnknown = false;
    CLIter it = children.begin();
    while(it != children.end()){
      assert((*it)->t == BTERM);
      CLIter e = it;
      it++;
      if ((*e)->bt == UNKNOWN_LITERAL){
	removedUnknown = true;
	(*e)->destroyTree();
	children.erase(e);
      }
    }
    if (!removedUnknown){
      this->destroyTree();
      return false;
    }
  }else if ( children.size()>3 ){
    // Note: all the children will be distinct, due to noDuplicateNodes()
    debout("cn4");
    this->destroyTree();
    return false;
  }//Children size == 2
  debout("cn4");

  // Sort the 2-tuple
  CLIter one,two;
  one = two = children.begin(); 
  two++;

  if( (*one)->bt > (*two)->bt ){
    BTerm temp = (*two)->bt; 
    (*two)->bt = (*one)->bt;
    (*one)->bt = temp;
  }
  return true;
}







SymNode *SymNode::copyTree(){
  SymNode *newRoot = new SymNode(this);

  CLIter iter, end;
  iter = children.begin();
  end = children.end();
  while(iter != end){
    newRoot->addChild((*iter)->copyTree());
    iter++;
  }
  return newRoot;
}
void SymNode::addChild( SymNode * nodePtr ){
  assert(nodePtr != 0);
  nodePtr->index = this->children.size();
  this->children.push_back(nodePtr);
}

void SymNode::addChildren(SymNode *other){
  CLIter iter, end;
  iter = other->children.begin();
  end  = other->children.end();
  while(iter != end){
    this->addChild( *iter );
    iter++;
  }
}
void SymNode::removeChild( int index ){
  //DOES NOTHING
  assert(0);
  assert( index >= 0 && (unsigned int)index < this->children.size() );
  int size = children.size();
  for(int i = 0; i < size; i++){
    
  }
  return;
}


void SymNode::removeParens(){
  // If child node has same type as parent, then we flatten the tree.
  CLIter cur = children.begin();
  while( cur != children.end() ){
    if((*cur)->t == t){
      debout("removeParens simplifies");
      this->addChildren(*cur);
      delete (*cur);
      cur = children.erase(cur);
    }else{
      cur++;
    }
  }
}

///////////////////////
// Post-condition (memory management):
//   Always return a pointer to this node.
//
SymNode* SymNode::toDNF(){
  debout("toDNF");
  CLIter iter,end,dup;
  if( t == BTERM || t == FTYPE || t == ETYPE || t == TTYPE){
    debout("p:isLiteral");
    return this;
  }
  //Ensure all of the clauses are in DNF
  iter = children.begin();
  end = children.end();

  ChildList *newList = new ChildList();
  while( iter != end ){
    debout("loop:p1");
    SymNode *returned;
    returned = (*iter)->toDNF();
    if( returned != NULL ){
      newList->push_back(*iter);
    }
    iter++;
  }
  debout("p2");
  this->children.swap( *newList );
  newList->clear();
  delete newList;

  deboutnn( "toDNF: children.size(): ");
  debout(children.size() );

  //All chilren are false if they have been pruned by DNF  
  if( children.empty() ){
    debout("p:shouldn\'t be happening");
    assert(0);
    this->destroyTree();
    return new SymNode(FTYPE,-1);
  }
  debout("p3");
  if( t == SOR ){
    debout("p3.1");
    
    //Then we are fine because all of the clauses are in DNF
    this->removeParens();
    //noDuplicateNodes(this->children);
    return this;
  } else if( t == SAND ){
    debout("p4");
    //We need to distribute "and" over the clauses that are in DNF already
    //It's exponential time!!!!!!!!!
    //TEMP MARK

    // BK: this never happens, because distribute has not line: return 0;
    if( this->distribute() ==  0 ){ 
      assert(0); //Shouldn't be happening
      debout("p6 this->distribute() == 0");
      //The children are false according to our restrictions and so this is false and therefor useless
      this->destroyTree();  // BK fixed possible leak
      delete this;
      return new SymNode(FTYPE,-1);
    }
    debout("p5");
    //noDuplicateNodes(this->children);
    return this;
  }else{
    //raise exception "not" is not implemented
    assert(0);
  }
  return NULL;
}


///////////////////
// Post-conditions (memory management): 
//    Every node listed in bterms is either affixed to this node or deleted.
//
// Distributes a conjunction of clauses.
//
int SymNode::distribute(){
  if(t != SAND){
    cerr << "distribute called on a non \"and\" node"<<endl;
    assert(t != SAND);
  }
  debout("d1");
  //Ensure that there are no ANDs as children
  removeParens();
  debout("d1.1");

  CLIter current;
  ChildList * bterms = new ChildList();
  current = children.begin();
  // Separate the complex terms from the BTerms;
  // children will contain only complex terms (i.e. ORs)
  while( current != children.end() ){
    if((*current)->t == BTERM){
      bterms->push_back(*current);
      current = children.erase(current);
    }else{
      current++;
    }
  }
 
  debout("d2");
  // Continue to call distr(Node,Node) until all of my children are 
  // distributed.  Recall that children contain ONLY complex nodes.
  while(children.size() > 1 ){
    debout( "d3.1");
    SymNode *result = distr( *(children.begin()), *(++children.begin()) );
    // BK 8/21/2009 cull to save run-time
    if (result->t != SAND){
      for(CLIter rc = result->children.begin(); rc != result->children.end(); ){
	if( (*rc)->cullNode2() ){
	  rc++;
	}else{
	  rc = result->children.erase(rc);
	}
      }
    }
    // END BK 8/21/2009
    children.pop_front();
    children.pop_front();
    debout("d3.2\nresult tree");
    deb(result->printTree(0));
    debout("d3.2.1");
    //MARK added this 
    assert(result!=0);
    if(result->t == SAND){
      // Result is an AND node, so add its children to the bterms list
      debout("result->t == SAND)");
      bterms->splice(bterms->begin(),result->children);
      noDuplicateNodes(*bterms);
      delete result;
      debout("leave result->t == SAND");
      debout("bterms size: "); 
      deboutnn(bterms->size());
    }else if( result->t == SOR ){
      // Result being an OR, means we added to the children list
      noDuplicateNodes(result->children);
      addChild(result);
    }else if(result->t == BTERM){
      bterms->push_back(result);
    }else if(result->t == FTYPE){
      delete result;
    }else{
      cerr << "a returned result was not an \"and\" or an \"or\" or a BTERM" <<endl;
      assert(0);
    }
  }

  debout("d4");
  if( children.size() == 0 && bterms->size() == 1 ){
    cerr << "not possible yet.1" <<endl;    
    assert(0);
    this->t = BTERM;
    this->bt = (*(bterms->begin()))->bt;
    delete bterms;  // BK fixed unreachable leak
    return 1;
  }else if( children.size() == 0 && bterms->size() > 1){
    
    //I stay an AND, but if we culled the ORs then I should be culled as well
    //That scenario isn't implemented yet
    children.swap(*bterms);
    delete bterms;
    return 1;
  }else if( children.size() == 0 && bterms->size() == 0){
    cerr << "not possible yet." <<endl;
    assert(0);
  }

  // Distribute the leaf bterms into the remaining complex term.
  debout("d5");
  if( children.size() == 1 ){
    this->t = SOR;
    this->bt = -1;
    SymNode *result = *(children.begin());
    children.pop_front();
    while (bterms->begin() != bterms->end()){
      SymNode* cur = *(bterms->begin());
      bterms->pop_front();
      result = distr(result, cur);
      // BK 8/21/2009 cull to save run-time
      for(CLIter rc = result->children.begin(); rc != result->children.end(); ){
        if((*rc)->cullNode2() ){
          rc++;
        }else{
          rc = result->children.erase(rc);
        }
      }
      noDuplicateNodes(result->children);
      // END BK 8/21/2009

      debout("distr returned a:");
      deb(result->printTree(0));
    }
    // copy the children of the resulting OR into the current children list
    this->addChildren( result );
    delete result;  
    delete bterms;   // BK fixed leak (8/7/2007)
    return 1;
  }else{
    cerr << "there should at least be one child and no bterms or no children and more than 1 bterm" <<endl;
    assert(0);
  }
    

  // BK: check the size of the bterms array pointer
  assert (bterms->size() == 0);


  delete bterms;   // BK fixed leak (8/7/2007)
  return 1;
}


///////////////////
//
// Distributes a conjunction (AND) over a pair of nodes.
//
// Post-conditions (memory management): 
//    Every newly created node is either in the return value or was deleted.
//    Every node (or copy of a node) that is not attached to the return value 
//    is deleted.
// Pre- & Post-condition (DNF sub-tree):
//    Return either an AND with at least two BTERM children, or an OR with 
//    every child being an AND with at least two children.  Meaning that
//    1 is always represented as (1 /\ 1) in an OR tree.
//
//
SymNode *distr(SymNode *n1, SymNode *n2){
  if (n1->t == ETYPE && n2->t == ETYPE)
    return n1;
  else if (n1->t == ETYPE)
    return n2;
  else if (n2->t == ETYPE)
    return n1;

  //Ensure that the first Predicate is always complex unless both are BTERMS
  if( n1->t == BTERM ){
    SymNode *temp = n2;
    n2 = n1;
    n1 = temp;
  }
  debout("distr is given:n1\n");
  deb(n1->printTree(0));
  debout("n2\n");
  deb(n2->printTree(0));
  debout("");
  // Because parent of these is an AND, and the code guarantees a "flat"
  // tree (i.e. no extra AND nodes).
  assert( (n1->t == SOR || n1->t == BTERM)  && (n2->t == SOR || n2->t == BTERM)  );

  deboutln("di1");
  if( n1->t == BTERM ){
    //Both are BTerms
    debout("di1.1");
    SymNode *newAnd = new SymNode(SAND, -1);
    newAnd->addChild(n1);
    newAnd->addChild(n2);
    // copied pointers to n1 or n2
    // CANNOT destroy them later
    return newAnd;
  }
  deboutln("di2");

  // CASE: n1 is OR and n2 is BTERM
  if(n2->t == BTERM){
    debout("di2.1");
    for(CLIter cur = n1->children.begin(); cur != n1->children.end(); cur++){
      SymNode *copyOfn2 = new SymNode(n2);
      if( (*cur)->t == BTERM){
	(*cur)->addChild(new SymNode(BTERM, (*cur)->bt));
	(*cur)->addChild( copyOfn2 );
	(*cur)->t = SAND;
	(*cur)->bt = -1;
      }else if((*cur)->t == SAND){
	(*cur)->addChild(copyOfn2);
      }else{
	assert((*cur)->t != SAND);
	copyOfn2->destroyNodeShallow(); // BK fixed possible leak
      }
    }

    // Distributed BTERM from n2 across all children of n1
    // made the appropriate number of copies of n2 and attached to n1
    // Changed node n1 to be type SAND
    // What happened to n2?  SHOULD DESTROY IT.
    n2->destroyNodeShallow(); // BK fixed leak (8/7/2007)
    n2 = NULL;                // BK fixed leak
    return n1;
    // returned an OR rooted subtree with AND children.
  }

  deboutln("di3");
  assert(n1->t == SOR);  
  assert(n2->t == SOR);  
  
  // BK 8/21/2009 culling while expanding
  // requires this case
  if (n2->children.size() == 0){
    n2->destroyNodeShallow();
    n2 = NULL;
    return n1;
  }

  // Distribute by making an AND node with children:
  //    1) a full copy of n1, and 2) a node of n2
  // Make the first distribution use the original of n1
  // And subsequent copies will use a copy of the original
  SymNode *newOr = new SymNode(SOR,-1);
  //Add the first "and" manually
  SymNode *newAnd = new SymNode(SAND, -1);
  deboutln(*(n2->children.begin()));
  newAnd->addChild( *(n2->children.begin()) ) ; // 831
  newAnd->addChild(n1);
  newOr->addChild(newAnd);
  deboutln("di4");

  for(CLIter cur = ++(n2->children.begin()); cur != n2->children.end(); cur++){
    //debout("di4.1");
    newAnd = new SymNode(SAND, -1);
    newAnd->addChild((*cur));
    newAnd->addChild(n1->copyTree()); // deep copy of n1
    newOr->addChild(newAnd);
  }
  n2->destroyNodeShallow(); // BK replace above with this (equivalent statement)
  n2 = NULL;
  debout("di5");
  // n1 is pointed to by newAnd, therefore, it should NOT be destroyed.


  //MARK changed to newOr
  for(CLIter cur = newOr->children.begin(); cur != newOr->children.end(); cur++){
    (*cur)->distribute();
  }
  newOr->removeParens(); // flatten extra ORs
  return newOr;
  // returned a flat OR sub-tree
}

int SymTree::get_next_pair(int &a, int &b){
  //Used to iterator over a culled tree
  //ready == 1 indicated just culled or able to give new node
  //ready == 0 means that there aren't any items left, the values for a,b are unchanged.
  if(ready == 1){
    //printf("SymTree::get_next_pair()\n");
    debout("g1");
    if( root->t == SAND ){
      debout("g1.1");
      assert(root->children.size() == 2);
      CLIter childI = root->children.begin();
      a = (*childI)->bt;
      childI++;
      b = (*childI)->bt;
      ready = 0;

      assert(a <= b); // BK: bug detection (8/7/2007)
      return 1;
    }else if(root->t == SOR){
      debout("g1.2");
      if((*curAnd)->t == BTERM){
	a = b = (*curAnd)->bt;
      }else{
	assert((*curAnd)->t == SAND);
	CLIter childI = (*curAnd)->children.begin();
	
	a = (*childI)->bt;
	childI++;
	b = (*childI)->bt;
      }
      curAnd++;

      if(curAnd == root->children.end()){
	ready = 0;
      }

      assert(a <= b); // BK: bug detection (8/7/2007)
      return 1;
    }else if(root->t == BTERM){
      debout("g1.3");
      a = b = root->bt;
      ready = 0;

      assert(a <= b); // BK: bug detection (8/7/2007)
      return 1;
    }else{
      return 0;
    }
  }else if(ready == 0){
    return 0;
  }else{//This was called after the tree was modified without culling
    cerr << "get_next_pair called without first calling cull. Status: "<<flush;
    cerr << ready << endl;
    assert(0);
    return ready;
  }
}



void SymTree::shuffle(){
  if(root->t == SOR && root->children.size() > 1){
    //printf("_UN_SHUFFLED TREE\n");
    //this->printTree();

    // copy the children over to old_order
    ChildList old_order = root->children;
    //for (CLIter d = root->children.begin(); d != root->children.end(); d++)
    //{
    //old_order.push_back(*d);
    //}
    root->children.clear();

    // uniformly at random (without replacement) 
    // add children from old_order back into root->children
    int num_children = old_order.size();
    for (int i= 0;  i < num_children;  i++)
      {
	double u = drand48(); // randomly selects the child to move
	double unit = 1.0  / (double) old_order.size(); // uniformly at random
	double prob = unit;
	for (CLIter c = old_order.begin(); c != old_order.end(); c++)
	  {
	    if ( u < prob )
	      {
		root->children.push_back(*c);
		old_order.erase(c);  // without replacement
		break;
	      }
	    prob += unit;
	  }  // end to select the child to move
	
      } // end for each child

    // print out nodes to see how well it worked
    //printf("SHUFFLED TREE\n");
    //this->printTree();

  } // end if root is SOR
}


void SymTree::rewind(){
  if(ready == 0 || ready == 1){
    ready = 1;
    if(root->t == SOR){
      curAnd = root->children.begin();
    }
    return;
  }
}


bool SymNode::equals(SymNode *other){
  // nodes must have same BTERMs or same number of children
  if(other->children.size() != children.size() 
      || !( t == other->t && bt == other->bt ) ){
    return false;
  }
  if(other->children.size() == 0 ){
    return true; 
  }

  CLIter tPtr = children.begin(); 
  CLIter oPtr = other->children.begin();
  while( tPtr != children.end() ){
    if( !(*tPtr)->equals(*oPtr) ){
      return false;
    }
    tPtr++;
    oPtr++;
  }
  return true;
}


// Compares all pairs of children, finding and removing duplicate BTERMs
void noDuplicateNodes( ChildList &l )
{
  //This is shallow!
  //ONLY SAFE WHEN DELETING BTERMS for now (note is prev to Aug 2009)

  // BK Aug 21, 2009: seems to be deep compare and deep delete
  ChildList dups;
  CLIter n1 = l.begin();
  CLIter n2 = ++(l.begin());
  CLIter end = l.end();
  while( n1 != end ){
    //assert((*n1)->t == BTERM);
    while( n2 != end ){
      //Comparison
      if( (*n1)->equals(*n2) ){
	dups.push_back(*n2);
      }
      n2++;
    }
    n1++;
    n2 = n1;
    n2++;
  }

  dups.sort();
  dups.unique();
  
  for(n2 = dups.begin(); n2 != dups.end(); n2++){
    for(n1 = l.begin(); n1 != l.end(); n1++){
      if(*n1 == *n2){
	(*n1)->destroyTree();
	l.erase(n1);
	break;
      }
    }
  }
  dups.clear();
}







int old_main( int args, char **argv ){
  /*
    Useful cull test cases:
    one literal T,F,BTerm
    one AND
    SOR in DNF
    SOR not in DNF
    SAND not in DNF
   */
  //To quiet compiler
  args = 1;
  char ** argv2 = argv;
  argv2 = argv2;

  //Test SymTree

  SymTree * st1 = new SymTree(1);
  st1 = and2(st1,new SymTree(2));
  st1 = or2(st1,3);
  st1 = or2(st1,4);

  SymTree * st2 = new SymTree(10);
  st2 = or2(st2,11);
  st2 = or2(st2,12);
  
  
  st1 = or2(st1,new SymTree(st2));
  cout << "AND and OR tests" << endl;
  st1->printTree();
  
  SymTree *st3;
  st3 = new SymTree(20);
  st3 = and2(st3,21);
  st3 = and2(st3,22);
  st3 = and2(st3,23);

  st3->printTree();
  delete st3;
 
  st1 = and2(st1,st2);

  st2 = 0;

  st1->printTree();
  //Should have been 20^21^22^(1\/2\/3\/4\/10\/11\/12)
  SymTree *tmp;
  cout << "TEST and2\n";
  tmp = and2(1,2);
  tmp->printTree();
  delete tmp;
  tmp = 0;

  cout << "TEST copy";
  SymTree *st1c = new SymTree(st1);
  st1 = or2(st1, and2(30,31) );
  st1c->printTree();
  
  delete st1c;
  st1c = 0;

  st1->printTree();
  
  delete st1;
  st1 = 0;

  cout << "TEST remove duplicate terms from SymNode"<<endl;

  SymNode * root = new SymNode(SAND, -1);
  root->addChild(new SymNode(BTERM, 1));
  noDuplicateNodes(root->children);
  cout <<"when only 1 node is present" <<endl;
  root->printTree(0);

  cout <<"when only 2 nodes are present, same." <<endl;
  root->addChild(new SymNode(BTERM, 1));
  noDuplicateNodes(root->children);
  root->printTree(0);

  cout <<"when 3 nodes are present, same." <<endl;
  root->addChild(new SymNode(BTERM, 1));
  root->addChild(new SymNode(BTERM, 1));
  root->addChild(new SymNode(BTERM, 1));
  noDuplicateNodes(root->children);
  root->printTree(0);
  

  cout <<"when 4 nodes are present, 3 same." <<endl;
  root->addChild(new SymNode(BTERM, 1));
  root->addChild(new SymNode(BTERM, 2));  
  root->addChild(new SymNode(BTERM, 1));
  root->addChild(new SymNode(BTERM, 1));
  noDuplicateNodes(root->children);
  root->printTree(0);

  root->destroyTree();
  root = 0;

  cout << "TEST toDNF "<< endl;

  SymTree *dnfTest = new SymTree(1);
  cout << "simple singleton\nbefore:"<< endl;
  dnfTest->printTree();
  dnfTest->toDNF();
  cout <<"after"<<endl;
  dnfTest->printTree();
  
  testIterator(dnfTest);

  cout << "simple 2-and\nbefore:"<< endl;
  dnfTest = and2( dnfTest,2);
  dnfTest->printTree();
  dnfTest->toDNF();
  cout <<"after"<<endl;
  dnfTest->printTree();
  
  testIterator(dnfTest);

  cout << "simple 3-and\nbefore:"<< endl;
  dnfTest = and2(dnfTest,3);
  dnfTest->printTree();
  dnfTest->toDNF();
  cout <<"after"<<endl;
  dnfTest->printTree();

  cout << "simple 4-and\nbefore:"<< endl;
  dnfTest = and2(dnfTest,4);
  dnfTest->printTree();
  dnfTest->toDNF();
  cout <<"after"<<endl;
  dnfTest->printTree();

  cout << "(5 or (and 1 2 3 4))\nbefore:"<< endl;  
  dnfTest = or2(dnfTest,5);
  dnfTest->printTree();
  dnfTest->toDNF();
  dnfTest->printTree();


  cout << "(or 5 6 (and 1 2 3 4))\nbefore:"<< endl;  
  dnfTest = or2(dnfTest,6);
  dnfTest->printTree();
  dnfTest->toDNF();
  dnfTest->printTree();
  
  cout << "(or 5 6 7 (and 1 2 3 4))\nbefore:"<< endl;  
  dnfTest = or2(dnfTest,7);
  dnfTest->printTree();
  dnfTest->toDNF();
  dnfTest->printTree();
  
  testIterator(dnfTest);

  cout << "(and 8 (or 5 6 ))\nbefore:"<< endl;  
  SymTree *dnfTestPrime = and2( new SymTree(8),or2(5,6));
  dnfTestPrime->printTree();
  dnfTestPrime->toDNF();
  dnfTestPrime->printTree();

  testIterator(dnfTestPrime);

  delete dnfTestPrime;
  dnfTestPrime = 0;

  cout << "(and 8 (or 5 6 7(and 1 2 3 4)))\nbefore:"<< endl;  
  dnfTest = and2(dnfTest,8);
  dnfTest->printTree();
  dnfTest->toDNF();
  dnfTest->printTree();
  

  cout << "(and 8 9 (or 5 6 7(and 1 2 3 4)))\nbefore:"<< endl;  
  //This below is a massive memory leak
  // BK comment: only leak in main (8/7/2007)
  SymTree *dnfTestSecond = and2(1,2)->sand(and2( 3,4))->sor(7)->sor(or2(5,6))->sand(and2(8,9));
  dnfTestSecond->printTree();
  cout << "after" <<endl;
  dnfTestSecond->toDNF();
  dnfTestSecond->printTree();


  cout << "(and 8 8 9 (or 5 6 7(and 1 2 3 4)))\nbefore:"<< endl;  
  dnfTestSecond = and2(dnfTestSecond,8);
  dnfTestSecond->printTree();
  cout << "after" <<endl;
  dnfTestSecond->toDNF();
  dnfTestSecond->printTree();
  
  cout << "\nCULLING Test:\nbefore" << endl;
  dnfTestSecond = or2(dnfTestSecond , and2( and2(13,12),and2(13,12)) );
  dnfTestSecond->printTree();
  dnfTestSecond->cull();
  cout << "after" << endl;
  dnfTestSecond->printTree();

  testIterator(dnfTest);

  delete dnfTest;
  dnfTest = NULL;

  delete dnfTestSecond;
  dnfTestSecond = 0;

  cout << "END TEST"<< endl;;


  // Success Case
  cout << "Success CASE" << endl;
  SymTree * tree = new SymTree();
  SymTree* children = new SymTree(4);
  SymTree* children2 = children->sand(5);
  delete children;
  children = children2;
  SymTree* tree2 = tree->sor(children);
  delete children;
  delete tree;
  tree=tree2;
  //delete children2;
  children = NULL;
  children = new SymTree(6);
  children2 = children->sand(7);
  delete children;
  children = children2;
  tree2 = tree->sor(children);
  delete children;
  delete tree;
  tree = tree2;
  children = new SymTree(6);
  children2 = children->sor(7);
  delete children;
  children = children2;
  tree2 = tree->sand(children);
  delete children;
  delete tree;
  tree = tree2;

  tree->printTree();
  cout << "toDNF" << endl;
  tree->toDNF();
  tree->printTree();
  cout << "cull" << endl;
  tree->cull();
  tree->printTree();
  cout << "END Success CASE" << endl << endl << endl;
  delete tree;



  SymTree *myTree = and2( and2(or2((BTerm)0,1),1),and2((BTerm)0,1));
  myTree->printTree();
		
  myTree->toDNF();
  myTree->printTree();
  myTree->cull();
  myTree->printTree();
  delete myTree;



  // Fail Case
  cout << "FAIL CASE:" << endl;
  SymTree* testCase = new SymTree();
  
  // Valid options (i.e. solutions)
  SymTree* myAnd = new SymTree();
  myAnd = and2((BTerm)0,5);
  testCase = or2(testCase, myAnd);
  myAnd = and2((BTerm)0,15);
  testCase = or2(testCase, myAnd);
  myAnd = and2((BTerm)0,7);
  testCase = or2(testCase, myAnd);
  myAnd = and2((BTerm)0,13);
  testCase = or2(testCase, myAnd);
  myAnd = and2(5,20);
  testCase = or2(testCase, myAnd);
  myAnd = and2(15,20);
  testCase = or2(testCase, myAnd);
  myAnd = and2(7,20);
  testCase = or2(testCase, myAnd);
  myAnd = and2(13,20);
  testCase = or2(testCase, myAnd);
  myAnd = and2(4,5);
  testCase = or2(testCase, myAnd);
  myAnd = and2(4,15);
  testCase = or2(testCase, myAnd);
  myAnd = and2(4,7);
  testCase = or2(testCase, myAnd);
  myAnd = and2(4,13);
  testCase = or2(testCase, myAnd);
  myAnd = and2(5,16);
  testCase = or2(testCase, myAnd);
  myAnd = and2(15,16);
  testCase = or2(testCase, myAnd);
  myAnd = and2(7,16);
  testCase = or2(testCase, myAnd);
  myAnd = and2(13,16);
  testCase = or2(testCase, myAnd);

  // All possible haplotype assignments for Father
//   SymTree* myAnd = new SymTree();
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)4,(BTerm)7);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)5,(BTerm)6);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)4,(BTerm)1);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)5,(BTerm)0);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)4,(BTerm)3);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)5,(BTerm)2);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)6,(BTerm)1);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)7,(BTerm)0);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)4,(BTerm)13);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)5,(BTerm)12);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)4,(BTerm)15);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)5,(BTerm)14);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)6,(BTerm)13);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)7,(BTerm)12);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)4,(BTerm)9);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)5,(BTerm)8);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)4,(BTerm)11);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)5,(BTerm)10);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)6,(BTerm)9);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)7,(BTerm)8);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)0,(BTerm)13);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)1,(BTerm)12);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)0,(BTerm)15);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)1,(BTerm)14);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)2,(BTerm)13);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)3,(BTerm)12);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)4,(BTerm)21);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)5,(BTerm)20);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)4,(BTerm)23);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)5,(BTerm)22);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)6,(BTerm)21);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)7,(BTerm)20);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)4,(BTerm)17);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)5,(BTerm)16);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)4,(BTerm)19);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)5,(BTerm)18);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)6,(BTerm)17);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)7,(BTerm)16);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)0,(BTerm)21);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)1,(BTerm)20);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)0,(BTerm)23);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)1,(BTerm)22);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)2,(BTerm)21);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)3,(BTerm)20);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)4,(BTerm)29);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)5,(BTerm)28);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)4,(BTerm)31);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)5,(BTerm)30);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)6,(BTerm)29);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)7,(BTerm)28);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)4,(BTerm)25);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)5,(BTerm)24);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)4,(BTerm)27);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)5,(BTerm)26);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)6,(BTerm)25);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)7,(BTerm)24);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)0,(BTerm)29);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)1,(BTerm)28);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)0,(BTerm)31);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)1,(BTerm)30);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)2,(BTerm)29);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)3,(BTerm)28);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)12,(BTerm)21);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)13,(BTerm)20);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)12,(BTerm)23);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)13,(BTerm)22);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)14,(BTerm)21);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)15,(BTerm)20);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)12,(BTerm)17);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)13,(BTerm)16);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)12,(BTerm)19);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)13,(BTerm)18);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)14,(BTerm)17);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)15,(BTerm)16);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)8,(BTerm)21);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)9,(BTerm)20);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)8,(BTerm)23);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)9,(BTerm)22);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)10,(BTerm)21);
//   testCase = or2(testCase, myAnd);
//   myAnd = and2((BTerm)11,(BTerm)20);
//   testCase = or2(testCase, myAnd);

  SymTree* testCase2 = new SymTree();
  SymTree* myOr = new SymTree();
  myOr = or2(or2(or2(5,15), 7), 13);
  testCase2 = and2(testCase2, myOr);
  myOr = or2(or2(or2(5,15), 7), 13);
  testCase2 = and2(testCase2, myOr);
  myOr = or2(or2(or2((BTerm) 0,20), 4), 16);
  testCase2 = and2(testCase2, myOr);
  myOr = or2(or2(or2(or2(or2(or2(or2((BTerm) 0, 21), 1), 20), 4), 17), 5), 16);
  testCase2 = and2(testCase2, myOr);
  myOr = or2(or2(or2(or2(or2(or2(or2(4, 15), 5), 14), 6), 13), 7), 12);
  testCase2 = and2(testCase2, myOr);
  myOr = or2(or2(or2(or2(or2(or2(or2(4, 15), 5), 14), 6), 13), 7), 12);
  testCase2 = and2(testCase2, myOr);
  cout << "  first step toDNF and cull" << endl;
  testCase2->toDNF();
  testCase2->cull();

  testCase = and2(testCase, testCase2);

  testCase->printTree();
  printf("BIG ToDNF\n");
  testCase->toDNF();
  printf("BIG Cull\n");
  testCase->cull();
  testCase->printTree();
  cout << "END FAIL CASE" << endl << endl << endl;

  delete testCase;


  return 0;
}

void testIterator(SymTree *tree){
  cout << "TEST iterator" << endl;

  tree->toDNF();
  tree->cull();

  cout <<"i1" << endl;

  int n = 3;
  int a = -100;
  int b = -100;
  while( n > 0 ){ 
    cout <<"[ ";
    while( tree->get_next_pair(a,b) == 1  ){
      cout << "(" << a << "," << b << ") ";
    }
    cout << "]" << endl;
    tree->rewind();
    n--;
  }
}
  //Below are little programs that I use to test C code
  /*
  int a = 1;
  int *b = &a;

  int &c = *b;
  a = 2;
  cout <<"\n";
  cout << c << endl;
  cout << "\n";

  */
  
  /*
  list<int> l;
  l.push_back(1);
  l.push_back(2);
  l.push_back(3);
  list<int>::iterator iter, dup;
  iter = l.begin();
  dup = iter;
  ++iter;
  l.erase(dup);
  cout << *(iter) << endl;
  */
  /*  int d = 1;
  int z= 10;
  int *e = &d;
  int *&f = e;
  cout << e << endl;
  e = &z;
  cout << f << endl;*/
  /*
  list<int> li;
  li.push_back(1);
  list<int>::iterator iter = li.begin();
  iter++;
  iter--;
  cout << *iter << endl;
  */
