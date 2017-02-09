#include "cramore.h"

/* cmd_vcf_delta_svm.cpp
 *
 * Copyright (C) 2015 Dongwon Lee and Hyun Min Kang
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

#include <cstdio>
#include <ctype.h>
#include <cstdlib>
#include <cstring>
#include <errno.h>

#include "libsvm_gkm.h"
#include "reference_sequence.h"
#include "bcf_ordered_reader.h"
#include "bcf_ordered_writer.h"

double calculate_lsgkm_score(svm_model* model, const char *seq)
{
    union svm_data x;
    double score;

    x.d = gkmkernel_new_object(seq, NULL, 0);

    svm_predict_values(model, x, &score);

    gkmkernel_delete_object(x.d);

    return score;
}

int32_t cmdVcfDeltaSVM(int32_t argc, char** argv) {
  std::string inVcf;
  std::string outVcf;
  std::string modelf;
  std::string infoKey("DELTA_SVM");
  std::string infoDesc;
  std::string modelList;
  std::string refFasta;
  int32_t window = 10;
  bool rawScore = false;
  int32_t verbosity = 10000;

  paramList pl;
  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input Options", NULL)
    LONG_STRING_PARAM("in-vcf",&inVcf,"Input VCF/BCF file")
    LONG_STRING_PARAM("model",&modelf,"lsgkm model file")
    LONG_STRING_PARAM("model-list",&modelList,"List of lsgkm model file, where the first column represents the INFO field key, and the second column represent the file path to the model file, and (optionally) third as description of each key")
    LONG_STRING_PARAM("ref",&refFasta,"Reference FASTA file")
    LONG_INT_PARAM("win",&window,"Window size of flanking sequence nearby the variant (on each side)")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out",&outVcf,"Output VCF/BCF file")
    LONG_STRING_PARAM("info",&infoKey,"INFO field key to add to the output VCF/BCF files")
    LONG_STRING_PARAM("desc",&infoDesc,"Description of the INFO key (to include in the BCF/VCF header)")    
    LONG_PARAM("raw-score",&rawScore,"Output raw lsgkm scores as [INFO.RAW] field")
    LONG_INT_PARAM("verbosity",&verbosity,"Verbosity of output")
  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();
  
  // sanity check of input arguments
  if ( inVcf.empty() || outVcf.empty() || refFasta.empty() ) {
    error("[E:%s:%d %s] --in-vcf, --out, --ref are required parameters",__FILE__,__LINE__,__FUNCTION__);
  }

  // load model list and open model files
  if ( modelf.empty() + modelList.empty() != 1 ) {
    error("[E:%s:%d %s] one and only of --model or --model-list options must be specified",__FILE__,__LINE__,__FUNCTION__);    
  }

  // read reference sequences
  ReferenceSequence ref(refFasta);

  // read model file, after loading, modelFH contains multiple models, infoKeys,Descs will contain INFO field values
  std::vector<std::string> modelFNs;
  //std::vector<svm_model*> modelFHs;
  std::vector<std::string> infoKeys;
  std::vector<std::string> infoDescs;
  if ( !modelList.empty() ) {  // if model list is not empty, we have multiple models
    htsFile* hp = hts_open(modelList.c_str(), "r"); // read the model list file
    if ( hp == NULL )
      error("[E:%s:%d %s] Cannot open file %s for reading",__FILE__,__LINE__,__FUNCTION__, modelList.c_str());
    kstring_t str = {0,0,0};
    int32_t lstr = 0;
    int32_t nfields = 0;
    int32_t* fields = NULL;
    // model list is assumed to have [INFO_KEY] [MODEL_FILE] [INFO_DESCRIPTION = INFO_KEY if empty]
    for( int32_t i=0; ( lstr = hts_getline(hp, KS_SEP_LINE, &str) ) >= 0; ++i ) {
      fields = ksplit(&str, 0, &nfields);
      if ( nfields < 2 )
	error("[E:%s:%d %s] Less than two columns observed in line %d of %s",__FILE__,__LINE__,__FUNCTION__, i+1, modelList.c_str());
      else if ( nfields > 3 )
	warning("[E:%s:%d %s] More than two columns observed in line %d of %s",__FILE__,__LINE__,__FUNCTION__, i+1, modelList.c_str());
      infoKeys.push_back(&str.s[fields[0]]);

      if ( nfields == 2 )
	infoDescs.push_back(&str.s[fields[0]]);
      else
	infoDescs.push_back(&str.s[fields[2]]);

      modelFNs.push_back(&str.s[fields[1]]);
      
      //svm_model* model = svm_load_model(&str.s[fields[1]]);  // load the model
      //if ( model == NULL )
      //	error("[E:%s:%d %s] Cannot load SVM model file %s for reading",__FILE__,__LINE__,__FUNCTION__, str.s[fields[1]]);
      //modelFHs.push_back(model);
    }
    hts_close(hp);
  }
  else { // model is not empty
    modelFNs.push_back(modelf);
    //svm_model* model = svm_load_model(modelf.c_str());
    //if ( model == NULL )
    //  error("[E:%s:%d %s] Cannot open file %s for reading",__FILE__,__LINE__,__FUNCTION__, modelf.c_str());
    //modelFHs.push_back(model);
    infoKeys.push_back(infoKey);

    if ( infoDesc.empty() ) 
      infoDescs.push_back(infoKey);
    else
      infoDescs.push_back(infoDesc);
  }

  // Read VCF file
  std::vector<GenomeInterval> intervals;
  BCFOrderedReader odr(inVcf.c_str(), intervals);
  BCFOrderedWriter odw(outVcf.c_str(), 0);
  bcf1_t* iv = bcf_init();

  odw.set_hdr(odr.hdr); // transfer header
  char buffer[65535];
  for(int32_t i=0; i < (int)infoKeys.size(); ++i) { // add INFO fields
    sprintf(buffer, "##INFO=<ID=%s,Number=A,Type=Float,Description=\"DeltaSVM for %s\">\n", infoKeys[i].c_str(), infoDescs[i].c_str());
    bcf_hdr_append(odw.hdr, buffer);
    if ( rawScore ) {
      sprintf(buffer, "##INFO=<ID=%s.RAW,Number=R,Type=Float,Description=\"Raw DeltaSVM for %s\">\n", infoKeys[i].c_str(), infoDescs[i].c_str());      
      bcf_hdr_append(odw.hdr, buffer);
    }
  }
  bcf_hdr_sync(odw.hdr);
  odw.write_hdr();

  if ( modelFNs.size() > 1 ) {
    std::vector< std::vector<std::string> > vseqs;
    std::vector<bcf1_t*> nvs;
    //std::vector< std::vector<float> > vscores;
    int32_t vcount;
    std::string alt;
    int32_t max_n_allele = 2;
    const char* rname;
    int32_t refLen;

    for(vcount=0; odr.read(iv); ++vcount) {
      nvs.push_back(bcf_dup(iv));
      // fetch allelic sequences and flanks
      bcf_unpack(iv, BCF_UN_STR);
      rname = bcf_hdr_id2name(odr.hdr, iv->rid); // rname is the sequence name
      
      if ( (vcount+1) % verbosity == 0 ) 
	notice("Reading %d variants from VCF at %s:%d. Beware that the memory size may increase dramatically as we load the whole VCF into memory", vcount+1, rname, iv->pos);

      vseqs.resize(vseqs.size()+1);
      vseqs.back().resize(iv->n_allele);
      refLen = iv->rlen > 2*window ? 2*window : iv->rlen;

      std::string& refSeq = vseqs.back()[0];
      
      ref.fetch_seq(rname, iv->pos+1-window, iv->pos+refLen+window, refSeq);
      for(int32_t i=0; i < (int32_t)refSeq.size(); ++i) {
	if ( refSeq[i] == 'N' ) refSeq[i] = 'A';
	switch( refSeq[i] ) {
	case 'A': case 'G': case 'C': case 'T':
	  break;
	case 'a': case 'g': case 'c': case 't':
	  refSeq[i] += ('A' - 'a');
	  break;
	default:
	  refSeq[i] = 'A';
	}
      }
      
      for(int32_t i=1; i < iv->n_allele; ++i) {
	std::string& altSeq = vseqs.back()[i];
	alt = iv->d.allele[i];
	for(int32_t j=0; j < (int32_t)alt.size(); ++j) {
	  if ( alt[j] == 'N' ) alt[j] = 'A';
	  switch( alt.at(j) ) {
	  case 'A': case 'G': case 'C': case 'T':
	    break;
	  case 'a': case 'g': case 'c': case 't':
	    alt.at(j) += ('A' - 'a');
	    break;
	  default:
	    alt.at(j) = 'A';
	  }
	}
	altSeq = refSeq.substr(0,window) + alt + refSeq.substr(window+refLen);
      }

      if ( max_n_allele < iv->n_allele )
	max_n_allele = iv->n_allele;
    }
    odr.close();

    notice("Finished loading %d variants", vcount+1);    

    for(int32_t i=0; i < (int32_t)modelFNs.size(); ++i) {
      notice("Processing SVM model file %s", modelFNs[i].c_str());
      svm_model* model = svm_load_model(modelFNs[i].c_str());

      std::vector<float> scores(max_n_allele);
      std::vector<float> deltas(max_n_allele-1);

      const char* infoKey = infoKeys[i].c_str();
      std::string infoKeyRawS = (infoKeys[i]+".RAW");
      const char* infoKeyRaw = infoKeyRawS.c_str();

      for(int32_t j=0; j < vcount; ++j) {
	for(int32_t k=0; k < nvs[j]->n_allele; ++k) {
	  scores[k] = (float)calculate_lsgkm_score(model, vseqs[j][k].c_str());
	  if ( k > 0 ) 
	    deltas[k-1] = scores[k] - scores[0];
	}

	bcf_update_info_float(odw.hdr, nvs[j], infoKey, deltas.data(), nvs[j]->n_allele-1);
	if ( rawScore )
	  bcf_update_info_float(odw.hdr, nvs[j], infoKeyRaw, scores.data(), nvs[j]->n_allele);

	if ( i == (int32_t)modelFNs.size()-1 ) {
	  odw.write(nvs[j]);
	  bcf_destroy(nvs[j]);
	}
      }

      svm_free_and_destroy_model(&model);
    }
  }
  else if ( modelFNs.empty() ) {
    error("[E:%s:%d %s] No SVM model file specified as argument --model or --model-list",__FILE__,__LINE__,__FUNCTION__);
  }
  else {
    svm_model* model = svm_load_model(modelFNs[0].c_str());  // load the model
    if ( model == NULL )
      error("[E:%s:%d %s] Cannot load SVM model file %s for reading",__FILE__,__LINE__,__FUNCTION__, modelFNs[0].c_str());
    
    // read each variant from VCF
    std::vector<std::string> seqs;
    std::vector<float> scores;
    std::vector<float> deltas;
    const char* rname;
    //int32_t cent1;
    std::string alt;
    int32_t vcount;
    int32_t refLen;
    
    for(vcount=0; odr.read(iv); ++vcount) {
      bcf1_t* nv = bcf_dup(iv);
      
      // fetch allelic sequences and flanks
      bcf_unpack(iv, BCF_UN_STR);
      rname = bcf_hdr_id2name(odr.hdr, iv->rid); // rname is the sequence name
      
      if ( (vcount+1) % verbosity == 0 ) 
	notice("Processing %d variants at %s:%d", vcount+1, rname, iv->pos);
      
      if ( (int32_t)seqs.size() < iv->n_allele ) {
	seqs.resize(iv->n_allele);
	scores.resize(iv->n_allele);
	deltas.resize(iv->n_allele-1);
      }
      // [pos+(rlen+1)/2-window,pos] . [pos+1, pos+rlen] . [pos+rlen+1, pos+(rlen+1)/2+window]
      // [0,window-(rlen+1)/2] . [window-(rlen+1)/2+1, window+(rlen-1)/2]. [window+(rlen+1)/2, 2*window]
      //cent1 = (iv->pos+(iv->rlen+1)/2);
      
      refLen = iv->rlen > 2*window ? 2*window : iv->rlen;
      
      // pos1 : 776546
      // iv->pos : 776545
      // refLen : 1
      // iv->pos+(refLen+1)/2-window+1 = 
      // iv->pos+(refLen+1)/2-window+1 = 776536
      // iv->pos+(refLen+1)/2+window+1 = 776556
      // refLen : 2
      // iv->pos+(refLen+1)/2-window+1 = 77653
      // iv->pos+(refLen+1)/2+window+1 = 776557    
      //ref.fetch_seq(rname, iv->pos+(refLen+1)/2-window+1, iv->pos+(refLen+1)/2+window+1, seqs[0]);
      ref.fetch_seq(rname, iv->pos+1-window, iv->pos+refLen+window, seqs[0]);
      for(int32_t i=0; i < (int32_t)seqs[0].size(); ++i) {
	if ( seqs[0][i] == 'N' ) seqs[0][i] = 'A';
	switch( seqs[0].at(i) ) {
	case 'A': case 'G': case 'C': case 'T':
	  break;
	case 'a': case 'g': case 'c': case 't':
	  seqs[0].at(i) += ('A' - 'a');
	  break;
	default:
	  seqs[0].at(i) = 'A';
	}
      }
      
      for(int32_t i=1; i < iv->n_allele; ++i) {
	alt = iv->d.allele[i];
	for(int32_t j=0; j < (int32_t)alt.size(); ++j) {
	  if ( alt[j] == 'N' ) alt[j] = 'A';
	  switch( alt.at(j) ) {
	  case 'A': case 'G': case 'C': case 'T':
	    break;
	  case 'a': case 'g': case 'c': case 't':
	    alt.at(j) += ('A' - 'a');
	    break;
	  default:
	    alt.at(j) = 'A';
	  }
	}
	seqs[i] = seqs[0].substr(0,window) + alt + seqs[0].substr(window+refLen);
	// window-refLen/2 = 10
	// window+(refLen+1)/2 = 11
	//if ( window > (refLen+1)/2 )
	//seqs[i] = seqs[0].substr(0,window-refLen/2) + alt + seqs[0].substr(window+(refLen+1)/2);
      }
      
      for(int32_t i=0; i < iv->n_allele; ++i) {
	scores[i] = (float)calculate_lsgkm_score(model, seqs[i].c_str());
      }
      
      // calculate delta-svm scores
      for(int32_t i=0; i < iv->n_allele-1; ++i) {
	deltas[i] = scores[i+1] - scores[0];
      }
      
      //for(int32_t i=0; i < iv->n_allele; ++i) {
      //notice("j=%d\ti=%d\tscore=%f",j,i,scores[i]);
      //}
      
      bcf_update_info_float(odw.hdr, nv, infoKeys[0].c_str(), deltas.data(), iv->n_allele-1);
      if ( rawScore )
	bcf_update_info_float(odw.hdr, nv, (infoKeys[0]+".RAW").c_str(), scores.data(), iv->n_allele);
      
      odw.write(nv);
      bcf_destroy(nv);
    }
    odr.close();
    
    notice("Finished processing %d variants", vcount+1);
    
    svm_free_and_destroy_model(&model);
  }

  odw.close();
  return 0;
}

