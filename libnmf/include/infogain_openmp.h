#ifndef UNIVIE_NMFLIB_INFOGAIN_OPENMP_H
#define UNIVIE_NMFLIB_INFOGAIN_OPENMP_H

#include "common.h"

void * infoGainThread(void * arg);

int infoGain_openmp(double  * data, int m, int n, double  * classColumn, int numDistClasses, idx_double  * indexedInfoGainSorted, idx_double  * indexedGainRatioSorted);

#endif