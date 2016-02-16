
#ifndef RECOMBPARSER
#define RECOMBPARSER




#include "genoparser.h"






class RecombParser : public GenoParser
{
 public:
  RecombParser(){
    fp = NULL;
  }


  void openRecomb(const char* recomb_file);
  void readRecomb(Pedigree& pedigree);

  void closeRecomb();


 private: 
  FILE* fp;
};



void RecombParser::openRecomb(const char* recomb_file)
{
  if (recomb_file != NULL)
    {
      fp = fopen (recomb_file, "r");
      if (fp == NULL)
	{
	  printf("Cannot open (r) %s.", recomb_file);
	  exit(-1);
	}
    }
}

void RecombParser::closeRecomb()
{
  if (fp != NULL)
    fclose(fp);
}


///////
// Reads whole recomb file.
void RecombParser::readRecomb(Pedigree& pedigree)
{


  // Each double appears on its own line.
  //


  char line [1024];
  int line_number = 1;
  int number_of_loci_read = 0;



  while (feof(fp) == 0)
     {
      
      // while not the end of the line
      char* notNull = fgets(line, 1024, fp);
      int length = strlen(line);
      char column[256];

      // Check for a blank line and break if there is...
      int check = 0;
      while (isspace(line[check]))
	    check++;
      if (check == length)
	break;


      if (notNull == NULL)
	break;

      //printf("READ: %s|\n", line);

      // read the recomb rate first
      double rec_rate = 1;
      if (sscanf(line, "%s", column) == 1)
	rec_rate = (double) atof(column);
      else
	printf("WARNING: Could not read the hap freq in %s\n", line);

      //printf("  rate %i] %f\n", line_number, rec_rate);
      pedigree.recombination_rates.push_back(rec_rate);


      number_of_loci_read++;



      // read through to an end of line character
      while (length > 0 && line[length-1] != '\n' && notNull != NULL)
	{
	  // read the next part of the buffer
	  //int length = strlen(line);
	  if (line[length-1] != '\n' && notNull != NULL){
	    notNull = fgets(line, 1024, fp);
	    length = strlen(line);
	  } else {
	    length = 0;
	  }
	  //printf("index %i, length %i, notNull %i \n", index, length, (notNull!=NULL));
	} // for each line



      line_number++;


      //if (number_of_loci_read == number_of_snps)
      //{
      //  // we finished reading the block
      //  break;
      //}
     
    } // while not end of file

     




  
}


#endif
