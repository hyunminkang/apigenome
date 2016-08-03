#ifndef UNIVIE_NMFLIB_INFOGAIN_SINGLEPREC_H
#define UNIVIE_NMFLIB_INFOGAIN_SINGLEPREC_H

#include "common.h"

void * infoGainThread(void * arg);

int infoGain_singleprec(float  * data, int m, int n, float  * classColumn, int numDistClasses, idx_float  * indexedInfoGainSorted, idx_float  * indexedGainRatioSorted);

#endif