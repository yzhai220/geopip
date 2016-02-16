

#ifndef GENOPARSER
#define GENOPARSER



#include <string>
#include <stdlib.h>





using namespace std;


class GenoParser
{
 public:
  GenoParser(int _smallest_allele, int _largest_allele){
    setAlleleBounds(_smallest_allele, _largest_allele);
  }

  void setAlleleBounds(int _smallest_allele, int _largest_allele);
  int readAllele(char* line, int length, int& index); // returns error code



  const static int ERR_OUTOFBUFFER; //OUTOFBOUNDS
  const static int ERR_INVALID_ALLELE;

  const static int UNKNOWN_ALLELE;


 protected:
  GenoParser(){}
  int smallest_allele;
  int largest_allele;

};

const int GenoParser::UNKNOWN_ALLELE = 0;

const int GenoParser::ERR_OUTOFBUFFER = -33;
const int GenoParser::ERR_INVALID_ALLELE = -22;






void GenoParser::setAlleleBounds(int _smallest_allele, int _largest_allele)
{
  if (_smallest_allele <= UNKNOWN_ALLELE){
    printf("FATAL ERROR in GenoParser(%i,%i): smallest allele <= %i\n", _smallest_allele, _largest_allele, UNKNOWN_ALLELE);
    exit(-1);
  }
  if (_smallest_allele > _largest_allele){
    printf("FATAL ERROR in GenoParser(%i,%i): smallest allele > largest allele\n", _smallest_allele, _largest_allele);
    exit(-1);
  }
  
  smallest_allele = _smallest_allele;
  largest_allele = _largest_allele;
}



int GenoParser::readAllele (char* line, int length, int& index)
{
  while (index < length && isspace(*(line + index)) ){
    index++;
  }
  //printf("read_allele(len %i, idx %i) line[idx] = %c\n", length, index, *(line + index));

  if (index >= length){ 
    //printf("Returning OUTOFBUFFER due to %i >= %i", index, length);
    return ERR_OUTOFBUFFER;
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
  //if (len > 1){
  //printf("double digit allele: |%s|\n", a.c_str());
  //return DOUBLEDIGIT;
  //}


  // check for arbitrary marker allele (i.e. within the given upper and lower bound)
  if (allele != UNKNOWN_ALLELE){
    if (allele < smallest_allele){
      return ERR_INVALID_ALLELE;
    }
    if (allele > largest_allele){
      return ERR_INVALID_ALLELE;
    }
  }

  while (index < length && isspace(*(line + index)) ){
    index++;
  }

  return allele;
}



#endif
