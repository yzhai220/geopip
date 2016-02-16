 #include "jni.h"
 #include <stdio.h>
 #include "pedi_io_Test.h"
 
 JNIEXPORT void JNICALL
 Java_pedi_io_Test_print(JNIEnv *env, jobject obj)
 {
     printf("Hello World from C!\n");
     return;
 }


JNIEXPORT jdoubleArray JNICALL Java_pedi_io_Test_iterate
  (JNIEnv *env, jobject obj, jdoubleArray arr)
{
  jsize size = (*env)->GetArrayLength(env, arr);

  printf("Array size: %i\n", size);

  double buf[size];
  double sum = 0;
  jint i;
  (*env)->GetDoubleArrayRegion(env, arr, 0, size, buf);
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


  double tmp[size]; /* make sure it is large enough! */
  int j;
  jdoubleArray iarr = (*env)->NewDoubleArray(env, size);
  if (iarr == NULL) {
    return NULL; /* out of memory error thrown */
  }
  for (j = 0; j < size; j++) {
    tmp[j] = buf[j] + 1;
  }
  (*env)->SetDoubleArrayRegion(env, iarr, 0, size, tmp);

  return iarr;
}
