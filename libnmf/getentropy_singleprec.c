// information gain C-Adaptation
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <sys/time.h>
#include <time.h>


#include "common.h"
#include "blaslapack.h"
#include "randnumber_singleprec.h"
#include "outputtiming.h"


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

//note numDistClasses is the number of rows of classvalueCount, in case of isClass = 0 => it is _not_ the number of unique values in classColumn
float  getEntropy_singleprec( float  * classColumn, int m, int isClass, float  * classvalueCount, int numDistClasses) {
  float  H = 0.0;		//Entropy
  int numInstances = m;		//number of instances
  int i, j, l; 			//loop index variables
  

  j = 1;
  classvalueCount[0] = classColumn[0];
  classvalueCount[numDistClasses] = 1.0;
  //copy list of distinct classes into first column of classvalueCount
  for (i = 1; i < m; ++i) {
    //test if classColumn[i] is already stored in classValueCount (up to row j)
    int unique = 1;
    for ( l = 0; l < j; ++ l) {
      if (classColumn[i] == classvalueCount[l]) {
	unique = 0;   
	classvalueCount[l + numDistClasses] += 1.0;
	break;
      }
    }
    if (unique) {
      classvalueCount[j] = classColumn[i];
      classvalueCount[j + numDistClasses] = 1.0;
      ++j;
    }
  }
  
  
  int actuallyUniqueValues = j;
  
  
  //for all unique class values
  for (i = 0; i < actuallyUniqueValues; ++i) {
    // incrementally calculate entropy
    // - sum(i=1 to k) pi log2 pi
    H = H - classvalueCount[i + numDistClasses]/numInstances * log2(classvalueCount[i + numDistClasses]/numInstances);
  }
  
  
  return H;
}