

#ifndef PEDIGREEPARSER
#define PEDIGREEPARSER



#include "genoparser.h"
#include "pedigree.h"





class PedigreeParser : public GenoParser
{
 public: 
  // Mendel algorithm needs sorted alleles to enumerate haplotypes
  // IBD HMM algorithm needs unsorted alleles to consider given haplotypes
  PedigreeParser(int _smallest_allele, int _largest_allele, bool _sortAlleles){
    setAlleleBounds(_smallest_allele, _largest_allele);
    sortAlleles = _sortAlleles;
  }

  
  void readPedigree(const char* in_file, Pedigree& pedigree);

 protected:
  bool sortAlleles;
};







void PedigreeParser::readPedigree(const char* in_file, Pedigree& pedigree)
{
  pedigree.setAlleleBounds(smallest_allele, largest_allele);


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

	int isFounder = false;
	if (faID == UNKNOWN_ALLELE && moID == UNKNOWN_ALLELE)
	  {
	    isFounder = true;
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
	fam = find(pedigree.families.begin(), pedigree.families.end(), newfam);
	if (fam == pedigree.families.end()){
	  // make family
	  pedigree.families.push_back(newfam);
	  fam = pedigree.families.end() - 1;
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
	int isTyped = false;
	
	
	char* notNull = fgets(line, 1024, fp);
	int length = strlen(line);
	int odd_allele = ERR_INVALID_ALLELE;

	//printf("line: %i\n", line_number);

	// in case it didn't read an end-of-line character
	while (length > 0 && notNull != NULL)
	  {
	    // for each SNP
	    int index = 0;
	    while (index  < length) // index++ occurs in readAllele
	      {
		int allele1 = ERR_INVALID_ALLELE;
		if (odd_allele != ERR_INVALID_ALLELE){
		  allele1 = odd_allele;
		  odd_allele = ERR_INVALID_ALLELE;
		} else {
		  //printf("1) readAllele(len %i, idx %i) line[idx] = |%c|\n", length, index, line[index]);
		  allele1 = readAllele(&line[0], length, index);
		}
		if (index >= length && allele1 != ERR_OUTOFBUFFER){
		  //printf("Odd allele: allele1 %i\n", allele1);
		  odd_allele = allele1;
		  break;
		}
		//printf("2) readAllele(len %i, idx %i) line[idx] = |%c|\n", length, index, line[index]);
		int allele2 = readAllele(&line[0], length, index);
		//printf("2) readAllele(len %i, idx %i) line[idx] = |%c|\n", length, index, line[index]);

		if (allele1 == INVALID_ALLELE || allele2 == INVALID_ALLELE){
		  printf(" INPUT ERROR: invalid allele (not btw %i and %i) read on line %i\n", smallest_allele, largest_allele, line_number);
		  exit(-1);
		}
		if ((allele1 == ERR_OUTOFBUFFER && allele2 != ERR_OUTOFBUFFER)
		     || (allele2 == ERR_OUTOFBUFFER && allele1 != ERR_OUTOFBUFFER)){
		  printf(" FATAL ERROR: odd number of alleles in fgets() on line %i\n", line_number);
		  printf("   allele1 %i   allele2 %i\n", allele1, allele2);
		  printf("   at locus # %i\n", ((int) (*person)->genotypes.size()+1) );
		  printf("   index %i, length %i\n", index, length);
		  printf("   odd_allele %i\n", odd_allele);
		  exit(-1);
		}

		if (allele1 != ERR_OUTOFBUFFER && allele2 != ERR_OUTOFBUFFER)
		  {
		    if (!isTyped && (allele1 != UNKNOWN_ALLELE || allele2 != UNKNOWN_ALLELE))
		      {
			isTyped = true;
		      }

		    //printf("  %i/%i", allele1, allele2); fflush(stdout);
		    Marker snp(allele1, allele2);
		    if (!sortAlleles)
		      {
			snp.allele1 = allele1;
			snp.allele2 = allele2;
		      }
		    
		    (*person)->genotypes.push_back(snp);
		    (*person)->inferred_genotypes.push_back(Marker(UNKNOWN_ALLELE, UNKNOWN_ALLELE));
		    (*person)->parental_source.push_back(UNKNOWN_ALLELE);

		    // add to list of distinct alleles
		    unsigned int index = (*person)->genotypes.size()-1;
		    while (pedigree.marker_alleles.size() <= index)
		      {
			// create new entries in marker_alleles
			pedigree.marker_alleles.push_back(DistinctAlleles());
		      }
		    if (allele1 != UNKNOWN_ALLELE)
		      pedigree.marker_alleles[index].insert(allele1);
		    if (allele2 != UNKNOWN_ALLELE)
		      pedigree.marker_alleles[index].insert(allele2);
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


#endif
