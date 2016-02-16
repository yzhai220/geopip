 #include "jni.h"
 #include <stdio.h>
 #include "pedi_hmm_Wrapper.h"

#include "pedigree.h"
#include "pedparser.h"
#include "recombparser.h"
#include "ibd_geno_alg.h"
 


Pedigree* pedigree = NULL;
IBDAlg ibd_alg;



 JNIEXPORT void JNICALL
 Java_pedi_hmm_Wrapper_print(JNIEnv *env, jobject obj)
 {
     printf("Hello World from C++!\n");
     return;
 }


JNIEXPORT jdoubleArray JNICALL Java_pedi_hmm_Wrapper_iterate
  (JNIEnv *env, jobject obj, jdoubleArray arr)
{
  //jsize size = (*env)->GetArrayLength(env, arr);
  jsize size = (jsize) env->GetArrayLength(arr);

  printf("Array size: %i\n", size);

  double buf[size];
  double sum = 0;
  jint i;
  //(*env)->GetDoubleArrayRegion(env, arr, 0, size, buf);
  env->GetDoubleArrayRegion(arr, 0, size, buf);
  for (i = 0; i < size; i++) {
    printf("%f ", buf[i]);
    sum += buf[i];
  }
  printf(" ");

  printf("sum: %f\n", sum);


  for (int k = 0;  k  < 5;  k++)
    {
      printf("blah...");
    }
  printf("\n");


  double tmp[size]; /* make sure it is large enough! */
  int j;
  //jdoubleArray iarr = (*env)->NewDoubleArray(env, size);
  jdoubleArray iarr = (jdoubleArray) env->NewDoubleArray(size);
  if (iarr == NULL) {
    return NULL; /* out of memory error thrown */
  }
  for (j = 0; j < size; j++) {
    tmp[j] = buf[j] + 1;
  }
  //(*env)->SetDoubleArrayRegion(env, iarr, 0, size, tmp);
  env->SetDoubleArrayRegion(iarr, 0, size, tmp);

  return iarr;
}



JNIEXPORT void JNICALL Java_pedi_hmm_Wrapper_loadfiles
  (JNIEnv *env, jobject obj, jstring jpedfile, jstring jrecombfile)
{
  // destruct pedigree if necessary
  if (pedigree != NULL)
    delete pedigree;
  pedigree = new Pedigree();


  // Convert Strings
  const char *recomb_file;
  const char *ped_file;
  //const jbyte *str;
  ped_file =  env->GetStringUTFChars(jpedfile, NULL);
  recomb_file = env->GetStringUTFChars(jrecombfile, NULL);
  if (ped_file == NULL || recomb_file == NULL) {
    return; /* OutOfMemoryError already thrown */
  }


  printf("READING input pedigree...\n");
  PedigreeParser ped_parser(1,2,0); // SNP DATA ONLY
  ped_parser.readPedigree(ped_file, *pedigree);
  pedigree->setFamilyPointers();
  printf("Traversal Order:\n");
  pedigree->setTraversalOrder(true);



  printf("READING input recombination rates...\n");
  RecombParser rec_parser;
  rec_parser.openRecomb(recomb_file);
  rec_parser.readRecomb(*pedigree);
  rec_parser.closeRecomb();


  // Release Strings
  env->ReleaseStringUTFChars(jpedfile, ped_file);
  env->ReleaseStringUTFChars(jrecombfile, recomb_file);

}


JNIEXPORT void JNICALL Java_pedi_hmm_Wrapper_setVerbose
  (JNIEnv *env, jobject obj, jboolean isVerbose)
{
  printf("Is verbose: %i\n", isVerbose);
  ibd_alg.setVerbose(isVerbose);
}


JNIEXPORT jdouble JNICALL Java_pedi_hmm_Wrapper_runHMM
  (JNIEnv *env, jobject obj)
{
  pedigree->runForwardBackwardBlockFamilyAlgorithm(ibd_alg, 1);

  vector<Family>::iterator fam;
  for (fam = pedigree->families.begin();  fam != pedigree->families.end();  fam++)
    {
      printf("Family %i\n", fam->family_id);
      printf("  Log Likelihood:  %f\n", fam->geno_log_prob_data);
    }

  fam = pedigree->families.begin();
  return  (jdouble) fam->geno_log_prob_data;
}



JNIEXPORT void JNICALL Java_pedi_hmm_Wrapper_setInhertianceIndicator
  (JNIEnv *env, jobject obj, 
   jint site, jint indiv_index, jint indiv_id, 
   jint parent, jdouble gf, jdouble gm)
{

  // set values only for first family, not all of them
  Family* family = &pedigree->families[0];

  if (family->members[indiv_index]->individual_id != indiv_id)
    {
      printf("ERROR: mismatch between indiv index and indiv id\n");
      fflush(stdout);
      jclass newExcCls;
      env->ThrowNew(newExcCls, "mismatch between indiv index and indiv id");
      return;
    }

  if (parent != 0 && parent != 1)
    {
      printf("ERROR: parent id must be zero or one\n");
      fflush(stdout);
      jclass newExcCls;
      env->ThrowNew(newExcCls, "ERROR: parent id must be zero or one\n");
      return;
    }

  Indicator indi(indiv_index, parent, gf, gm);


  int num_sites = family->members[0]->genotypes.size();
  if (family->indicator_weights.size() < num_sites)
    {
      family->indicator_weights.resize(num_sites);
    }

  family->indicator_weights[site].push_back(indi);
}



JNIEXPORT jdoubleArray JNICALL Java_pedi_hmm_Wrapper_getInheritanceIndicator
  (JNIEnv *env, jobject obj, jint site, jint indiv_index, jint indiv_id, 
   jint parent)
{
  // get values only for first family, not all of them
  Family* family = &pedigree->families[0];

  if (family->members[indiv_index]->individual_id != indiv_id)
    {
      printf("ERROR: mismatch between indiv index and indiv id\n");
      fflush(stdout);
      jclass newExcCls;
      env->ThrowNew(newExcCls, "mismatch between indiv index and indiv id");
      return NULL;
    }

  if (parent != 0 && parent != 1)
    {
      printf("ERROR: parent id must be zero or one\n");
      fflush(stdout);
      jclass newExcCls;
      env->ThrowNew(newExcCls, "ERROR: parent id must be zero or one\n");
      return NULL;
    }


   int size = 2;
  double grandparents[size]; /* make sure it is large enough! */

  grandparents[0] = family->marginal_indicators[site][4*indiv_index + 2*parent + 0];
  grandparents[1] = family->marginal_indicators[site][4*indiv_index + 2*parent + 1];

  // Make java array
  jdoubleArray return_array = (jdoubleArray) env->NewDoubleArray(size);
  if (return_array == NULL) {
    return NULL; /* out of memory error thrown */
  }
  //env->SetDoubleArrayRegion(return_array, 0, size, grandparents);

  return return_array;

}
