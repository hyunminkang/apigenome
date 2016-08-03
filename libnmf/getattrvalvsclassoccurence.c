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
void getattrValVsClassOccurence(double * attrCol, int m, double * classCol, int numDistClasses, double * classvalueCount, double * uniqueAttrVals, int * numUniqueAttrVals, int * classCountForAttrIndex) {
  int i, j; 			//loop variables
  
  //unique number of attribute values for this attribute
  *numUniqueAttrVals = 1;
  uniqueAttrVals[0] = attrCol[0];
  for (i = 1; i < m; ++i) {
    //check if attrCol[i] is not already stored in uniqueAttrVals
    int unique = 1;
    for (j = 0; j < *numUniqueAttrVals; ++j) {
      if (attrCol[i] == uniqueAttrVals[j]) {
	unique = 0;
	break;
      }
    }
    if (unique) {
      uniqueAttrVals[*numUniqueAttrVals] = attrCol[i];
      ++(*numUniqueAttrVals);
    }
  }
  

  for (i = 0; i < m*numDistClasses; ++i)
    classCountForAttrIndex[i] = 0;
  
  //for all instances do
  for (i = 0; i < m; ++i) {
    //find corresponding index of attrCol[i] in uniqueAttrVals
    int uniqi = 0;
    for (j = 0; j < *numUniqueAttrVals; ++j) {
      if (attrCol[i] == uniqueAttrVals[j]) {
	uniqi = j;
	break;
      }
    }
    
    //find corresponding classCol index
    int classIndex = 0;
    for (j = 0; j < numDistClasses; ++j) {
      if (classCol[i] == classvalueCount[j]) {
	classIndex = j;
	break;
      }
    }
    
    ++classCountForAttrIndex[uniqi + m * classIndex];
  }
  
}