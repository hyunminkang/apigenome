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
#include "randnumber.h"
#include "outputtiming.h"
#include "calculatenorm.h"
#include "calculatemaxchange.h"
#include "getattrvalvsclassoccurence.h"


/**
  * InfoGain calculates the InformationGain and GainRatio of the matrix "data" according to the classification given by classColumn
  * 
  * Parameter
  * data		double*		IN, matrix to calculate IG and GR for, instance times attribute matrix, attributes stored column wise
  * m			int		IN, first dimension of matrix data (rows)
  * n			int		IN, second dimension of matrix data (columns)
  * classColumn 	int*		IN, array of class attributes for every entity (row) in data
  * numDistClasses	int		IN, number of distinct classes for the instances in matrix "data"
  * indexInfoGain	int*		OUT, array of index of attributes sorted descendingly wrt information gain
  * infoGainSorted	double*		OUT, array of the according information gains of these attributes
  * indexGainRatio	int*		OUT, arry of index of attributes sorted descendingly wrt gain ratio
  * gainRatioSorted	double*		OUT, array of the according gain ratio of these attributes
  *
  */

//note numDistClasses is the number of rows of classvalueCount, in case of isClass = 0 => it is _not_ the number of unique values in classColumn
double getInfoGainForAttribute(double * attrCol, double * classCol, int numDistClasses, double * classvalueCount, int numInstances, double H) {
  double * uniqueAttrVals = (double*) malloc(sizeof(double) * numInstances);
  int * classCountForAttrIndex = (int*) malloc(sizeof(int) * numInstances * numDistClasses);
  int numUniqueAttrVals = 0;
  getattrValVsClassOccurence(attrCol, numInstances, classCol, numDistClasses, classvalueCount, uniqueAttrVals, &numUniqueAttrVals, classCountForAttrIndex);

  //allocate memor for CondEntropy array
  double * condEntropy = malloc(sizeof(double) * numUniqueAttrVals);
  
  
  
  
  //for all unique attribute values
  int i, j;				//loop variables
  int s = 0;				//total occurences of a attribute Value in attrCol
  for (i = 0; i < numUniqueAttrVals; ++i) { 
    //initalize condEntropy[i] to zero
    condEntropy[i] = 0.0;
    
    //attribute has value uniqueAttrVals(i) for s instances
    s = 0;
    for(j = 0; j < numDistClasses; ++j)
      s += classCountForAttrIndex[i + j * numInstances];
    
    
    //for all unique class values
    for(j = 0; j < numDistClasses; ++j) {
      //prevent division by zero
      if (classCountForAttrIndex[i + j * numInstances] != 0)
	condEntropy[i] = condEntropy[i] - classCountForAttrIndex[i + j * numInstances] / (double) s * log2( classCountForAttrIndex[i + j * numInstances] / (double) s); 
    }
    
    //compute conditional entropy for attribute value uniqueAttrVals[i]
    condEntropy[i] = (double) s / numInstances * condEntropy[i];
  }
  
  //compute the overall information gain for current attribute
  double infoGain = H;
  for (i = 0; i < numUniqueAttrVals; ++i) {
    infoGain -= condEntropy[i];
  }
  
  
  //free allocated memory
  free(condEntropy);
  free(classCountForAttrIndex);
  free(uniqueAttrVals);
  
  
  return infoGain;
}
