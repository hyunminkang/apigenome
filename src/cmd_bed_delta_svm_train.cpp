#include "cramore.h"

/* cmd_vcf_delta_svm.cpp
 *
 * Copyright (C) 2016 Dongwon Lee and Hyun Min Kang
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include "libsvm_gkm.h"
#include "libsvm.h"
#include "reference_sequence.h"

int32_t cmdBedDeltaSVMTrain(int32_t argc, char** argv) {
  int32_t kernelFunction = EST_TRUNC_PW;
  int32_t wordLength = 11;
  int32_t numInfoColumn = 7;
  int32_t maxMismatch = 3;
  double  gammaRBF = 1.0;
  int32_t cwMax = 50;
  double  cwHalfLife = 50.0;
  double  svmC = 1.0;
  double  epsilon = 1e-3;
  double  svmW = 1.0;
  double  cacheMB = 100;
  bool    shrinkHeuristic = false;
  int32_t nCrossVal = 0;
  int32_t idxCrossVal = 0;
  int32_t seed = 0;
  int32_t verbose = 1000;
  int32_t numThreads = 1;
  
  std::string posBed;
  std::string negBed;
  std::string refFasta;
  std::string outPrefix;

  paramList pl;
  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Required parameters", NULL)
    LONG_STRING_PARAM("pos",&posBed,"BED file representing the positive labels")
    LONG_STRING_PARAM("neg",&negBed,"BED file representing the negative labels")    
    LONG_STRING_PARAM("out",&outPrefix,"Output prefix")
    LONG_STRING_PARAM("ref",&refFasta,"FASTA format reference genome")    

    LONG_PARAM_GROUP("General parameters for LSGKM", NULL)
    LONG_INT_PARAM("kernel", &kernelFunction, "0: Gapped k-mer, 1: EST-FULL, 2: GKM, 3: GKM-RBF, 4: GKM-CenterWeighted, 5: GKM-CenterWeight-RBF")
    LONG_INT_PARAM("length",&wordLength, "Word length, between 3 and 12")
    LONG_INT_PARAM("k",&numInfoColumn, "Number of informative column. Should be less than length")
    LONG_INT_PARAM("max-dist",&maxMismatch, "Maximum mismatch to consider")
    LONG_INT_PARAM("threads",&numThreads,"Number of threads")

    LONG_PARAM_GROUP("Parameters for RBF kernels (3,5)", NULL)    
    LONG_DOUBLE_PARAM("gamma",&gammaRBF, "Gamma for RBF kernels")

    LONG_PARAM_GROUP("Parameters for center-weighted (WGKM) kernels (4,5 only)", NULL)        
    LONG_INT_PARAM("cw-max", &cwMax, "Maximum weight for WGKM kernels")
    LONG_DOUBLE_PARAM("cw-half", &cwHalfLife, "Half-life parameter of exponential decay")

    LONG_PARAM_GROUP("Standard libSVM parameters", NULL)        
    LONG_DOUBLE_PARAM("svm-c", &svmC, "Regularization (C) parameter in SVM")
    LONG_DOUBLE_PARAM("svm-w", &svmW, "Weight (W) parameter in SVM")    
    LONG_DOUBLE_PARAM("epsilon", &epsilon, "Precision parameters in SVM")
    LONG_DOUBLE_PARAM("cache-mb", &cacheMB, "Cache size in MB")
    LONG_PARAM("shrink", &shrinkHeuristic, "Turn on shrinking heutistic option")
    LONG_INT_PARAM("cv-num", &nCrossVal, "Number of cross-validation to perform (0: OFF)")
    LONG_INT_PARAM("cv-idx", &idxCrossVal, "1-based index of cross-validation run (when running in parallel)")

    LONG_PARAM_GROUP("Standard libSVM parameters", NULL)            
    LONG_INT_PARAM("seed", &seed, "Seed for random number generator")
    LONG_INT_PARAM("verbose",&verbose, "Verbosity interval")
   END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  // sanity check of input arguments
  if ( posBed.empty() || negBed.empty() || refFasta.empty() || outPrefix.empty() ) 
    error("[E:%s:%d %s] --pos --neg, --out, --ref are required parameters",__FILE__,__LINE__,__FUNCTION__);

  char abuf[65536];

  std::string modelFile;
  std::string predFile;
  if ( nCrossVal > 0 ) {
    if ( idxCrossVal > 0 ) 
      sprintf(abuf, "%s.cvpred.%d.txt", outPrefix.c_str(), idxCrossVal);
    else
      sprintf(abuf, "%s.cvpred.txt", outPrefix.c_str());
    predFile = abuf;
  }

  sprintf(abuf, "%s.model.txt", outPrefix.c_str());
  modelFile = abuf;

  struct svm_parameter param;
  param.svm_type = C_SVC;
  param.kernel_type = kernelFunction;
  param.L = wordLength;
  param.k = numInfoColumn;
  param.d = maxMismatch;
  param.M = cwMax;
  param.H = cwHalfLife;
  param.gamma = gammaRBF;
  param.cache_size = cacheMB;
  param.C = svmC;
  param.eps = epsilon;
  param.shrinking = shrinkHeuristic ? 1 : 0;
  param.nr_weight = 0;
  param.weight_label = new int[1];
  param.weight = new double[1];
  param.p = 0.1;
  param.probability = 0;
  param.nu = 0.5;

  gkmkernel_init(&param);
  gkmkernel_set_num_threads(numThreads);

  // read reference sequences
  ReferenceSequence ref(refFasta);

  std::vector<std::string> rnames;
  std::vector<int32_t> beg0s; // 30,32 .. 31,32 
  std::vector<int32_t> end1s;
  std::vector<bool> isPos;

  // load the BED files
  htsFile* hp = hts_open(posBed.c_str(), "r"); // read the model list file
  if ( hp == NULL )
    error("[E:%s:%d %s] Cannot open file %s for reading",__FILE__,__LINE__,__FUNCTION__, posBed.c_str());
  kstring_t str = {0,0,0};
  int32_t lstr = 0;
  int32_t nfields = 0;
  int32_t* fields = NULL;
  // model list is assumed to have [INFO_KEY] [MODEL_FILE] [INFO_DESCRIPTION = INFO_KEY if empty]
  for( int32_t i=0; ( lstr = hts_getline(hp, KS_SEP_LINE, &str) ) >= 0; ++i ) {
    fields = ksplit(&str, 0, &nfields);
    if ( nfields < 3 )
	error("[E:%s:%d %s] Less than three columns observed in line %d of %s",__FILE__,__LINE__,__FUNCTION__, i+1, posBed.c_str());
    
    rnames.push_back(&str.s[fields[0]]);
    beg0s.push_back(atoi(&str.s[fields[1]]));
    end1s.push_back(atoi(&str.s[fields[2]]));
    isPos.push_back(true);
  }
  hts_close(hp);

  hp = hts_open(negBed.c_str(), "r"); // read the model list file
  if ( hp == NULL )
    error("[E:%s:%d %s] Cannot open file %s for reading",__FILE__,__LINE__,__FUNCTION__, posBed.c_str());
  for( int32_t i=0; ( lstr = hts_getline(hp, KS_SEP_LINE, &str) ) >= 0; ++i ) {
    fields = ksplit(&str, 0, &nfields);
    if ( nfields < 3 )
	error("[E:%s:%d %s] Less than three columns observed in line %d of %s",__FILE__,__LINE__,__FUNCTION__, i+1, posBed.c_str());

    rnames.push_back(&str.s[fields[0]]);
    beg0s.push_back(atoi(&str.s[fields[1]]));
    end1s.push_back(atoi(&str.s[fields[2]]));
    isPos.push_back(false);
  }
  hts_close(hp);

  notice("Finished loading %u positive or negative labels", isPos.size());

  struct svm_problem prob;
  prob.l = (int32_t)isPos.size();
  prob.y = new double[prob.l];
  prob.x = new union svm_data[prob.l];

  std::string seq;

  for(int32_t i=0; i < prob.l; ++i) {
    if ( i % 10000 == 0 ) notice("Fetching %d sequences.. %s:%d-%d", i, rnames[i].c_str(), beg0s[i]+1, end1s[i]);
    if ( end1s[i] - beg0s[i] > MAX_SEQ_LENGTH ) {
      notice("%d-th label %s:%d-%d is too long (%d>%d bp). Using only %d bp in the the center", i+1, rnames[i].c_str(), beg0s[i]+1, end1s[i], end1s[i]-beg0s[i], MAX_SEQ_LENGTH, MAX_SEQ_LENGTH);
      int32_t mid = (end1s[i] + beg0s[i]) / 2;
      end1s[i] = mid + MAX_SEQ_LENGTH/2;
      beg0s[i] = mid + MAX_SEQ_LENGTH/2;
    }
    
    ref.fetch_seq(rnames[i], beg0s[i]+1, end1s[i], seq);
    //notice("i = %d seq = %s", i, seq.c_str());    
    prob.y[i] = isPos[i] ? 1 : -1;
    prob.x[i].d = gkmkernel_new_object(seq.c_str(), NULL, i);
    //notice("kernel created");    
  }

  const char* errorMsg = svm_check_parameter(&prob, &param);
  if ( errorMsg ) {
    error(errorMsg);
    exit(1);
  }

  notice("Starting training...");  
 
 if ( nCrossVal > 0 ) {
    srand(seed);
    double* target = new double[prob.l];
    svm_cross_validation(&prob, &param, nCrossVal, idxCrossVal, target, predFile.c_str());
    delete[] target;
  }
  else {
    svm_model* model = svm_train(&prob, &param);
    if ( svm_save_model(modelFile.c_str(), model) ) {
      error("Cannot save model to file %s",modelFile.c_str());
      exit(1);
    }
    svm_free_and_destroy_model(&model);
  }

  for(int32_t i=0; i < prob.l; ++i)
    gkmkernel_delete_object(prob.x[i].d);

  svm_destroy_param(&param);
  delete[] prob.y;
  delete[] prob.x;

  return 0;
}

