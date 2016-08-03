// information gain C-Adaptation
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <sys/time.h>
#include <time.h>
#include <omp.h>


#include "common.h"
#include "cpmidx_float.h"
#include "blaslapack.h"
#include "randnumber_singleprec.h"
#include "outputtiming.h"
#include "getentropy_singleprec.h"
#include "getinfogainforattribute_singleprec.h"
#include "infogain_openmp_singleprec.h"


/**
  * InfoGain calculates the InformationGain and GainRatio of the matrix "data" according to the classification given by classColumn
  * 
  * Parameter
  * data		float *		IN, matrix to calculate IG and GR for, instance times attribute matrix, attributes stored column wise
  * m			int		IN, first dimension of matrix data (rows)
  * n			int		IN, second dimension of matrix data (columns)
  * classColumn 	int*		IN, array of class attributes for every entity (row) in data
  * numDistClasses	int		IN, number of distinct classes for the instances in matrix "data"
  * indexInfoGain	int*		OUT, array of index of attributes sorted descendingly wrt information gain
  * infoGainSorted	float *		OUT, array of the according information gains of these attributes
  * indexGainRatio	int*		OUT, arry of index of attributes sorted descendingly wrt gain ratio
  * gainRatioSorted	float *		OUT, array of the according gain ratio of these attributes
  *
  */

int infoGain_openmp(float  * data, int m, int n, float  * classColumn, int numDistClasses, idx_float  * indexedInfoGainSorted, idx_float  * indexedGainRatioSorted) {
 
  
  //check number of supported cores
  int numProcs = omp_get_num_procs();
  int i;			//index variable used in for-loops
  int numInstances = m;		//number of instances in matrix data
  int numAttributes = n;	//number of attributes in matrix data
  float  H = 0.0;		//entropy
  float  * classvalueCount = (float *) malloc(sizeof(float ) * numDistClasses * 2);
  float  * attributevalueCount = (float *) malloc(sizeof(float ) * m * 2 * numProcs);
  
  // get Entropy of the classification 
  H = getEntropy_singleprec(classColumn, m, 1, classvalueCount, numDistClasses);


  
  float  * timing = (float *) malloc(sizeof(float ) * numAttributes * 1);	//timing in seconds
  float  splitinfo = 0.0;
  struct timeval start, end;
  
  
  
#pragma omp parallel num_threads(numProcs)
  {
  
#pragma omp for private(i, start, end, splitinfo)
  for (i = 0; i < numAttributes; ++i) {	
    
    int whoami = omp_get_thread_num();
    
    gettimeofday(&start, 0);
    
    indexedInfoGainSorted[i].val = getInfoGainForAttribute_singleprec(data + i * m, classColumn, numDistClasses, classvalueCount, numInstances, H);
//TODELETE
//printf("indexedInfoGainSorted[%d].val = %f\n", i, indexedInfoGainSorted[i].val);
    indexedInfoGainSorted[i].idx = i;
    //getEntropy(attributeColumn, noClassification, null) --- null since classvalueCount is not needed
    //TODO allocate second variable "attributevalueCount" of size m times 2
    splitinfo = getEntropy_singleprec(data + i * m, m, 0, attributevalueCount + whoami * 2 * numProcs, m);
    if (splitinfo != 0.0) {
      indexedGainRatioSorted[i].val = indexedInfoGainSorted[i].val / splitinfo;
    }
    else {
      indexedGainRatioSorted[i].val = 0.0;
    }
    indexedGainRatioSorted[i].idx = i;
    gettimeofday(&end, 0);
    // save timing in seconds
    timing[i] = (float ) ( ( 1.0 * end.tv_sec - start.tv_sec) + 1E-6 * (1.0 * (end.tv_usec - start.tv_usec)) );

  }
  
  }
  
  //TODO sort infoGain and GainRatio and store new indices of those values in descending order
  qsort((void*) indexedInfoGainSorted, numAttributes, sizeof(idx_float ), cmpidx_float );
  qsort((void*) indexedGainRatioSorted, numAttributes, sizeof(idx_float ), cmpidx_float );
  
  
  //free allocated memory
  free(classvalueCount);
  free(attributevalueCount);
  free(timing);
  
  // return 0 if everything worked fine, else return error code
  return 0;
}
