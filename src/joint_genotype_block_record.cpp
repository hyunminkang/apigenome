/* The MIT License
   Copyright (c) 2014 Adrian Tan <atks@umich.edu>
   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:
   The above copyright notice and this permission notice shall be included in
   all copies or substantial portions of the Software.
   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
   THE SOFTWARE.
*/

#include "joint_genotype_block_record.h"

/**
 * Constructor.
 * @v - VCF record.
 */
JointGenotypeBlockRecord::JointGenotypeBlockRecord(bcf_hdr_t *h, bcf1_t *v, int32_t vtype, int32_t nsamples, bool _printTmpInfo)
{
    clear();

    printTmpInfo = _printTmpInfo;

    //est = &globEst; //new Estimator();  
    
    this->h = h;
    //this->v = v;
    this->rid = bcf_get_rid(v);
    this->pos1 = bcf_get_pos1(v);
    this->vtype = vtype;
    this->nsamples = nsamples;
    
    this->alleles.l = this->alleles.m = 0; this->alleles.s = NULL;
    char** tmp_alleles = bcf_get_allele(v);
    for (size_t i=0; i< bcf_get_n_allele(v); ++i) {
      if (i) kputc(',', &this->alleles);
      kputs(tmp_alleles[i], &this->alleles);
      v_alleles.push_back(tmp_alleles[i]);
    }

    n_filter = 0;

    if (vtype==VT_SNP && v_alleles.size()==2)
    {
      //rid = bcf_get_rid(v);
      beg1 = bcf_get_pos1(v);
      end1 = beg1;

      if (bcf_has_filter(h, v, const_cast<char*>("overlap_snp"))==1) 
	n_filter |= FILTER_MASK_OVERLAP_SNP;
      
      if (bcf_has_filter(h, v, const_cast<char*>("overlap_indel"))==1) 
	n_filter |= FILTER_MASK_OVERLAP_INDEL;

      if (bcf_has_filter(h, v, const_cast<char*>("overlap_vntr"))==1) 
	n_filter |= FILTER_MASK_OVERLAP_VNTR;      
    }
    else if (vtype==VT_INDEL && v_alleles.size()==2)
    {
      //rid = bcf_get_rid(v);
      dlen = strlen(tmp_alleles[1])-strlen(tmp_alleles[0]);
      len = abs(dlen);
      
      int32_t *flanks_pos1 = NULL;
      int32_t n = 0;
      
      if (bcf_get_info_int32(h, v, "FLANKS", &flanks_pos1, &n)>0) {
	this->beg1 = flanks_pos1[0];
	this->end1 = flanks_pos1[1];
	free(flanks_pos1);
      }
      else {
	this->beg1 = bcf_get_pos1(v) - 3;
	this->end1 = bcf_get_end_pos1(v) + 3;
      }
      
      if (dlen>0) {
	indel.append(&tmp_alleles[1][1]);
      }
      else {
	indel.append(&tmp_alleles[0][1]);
      }

      if (bcf_has_filter(h, v, const_cast<char*>("overlap_snp"))==1) 
	n_filter |= FILTER_MASK_OVERLAP_SNP;
      
      if (bcf_has_filter(h, v, const_cast<char*>("overlap_indel"))==1) 
	n_filter |= FILTER_MASK_OVERLAP_INDEL;

      if (bcf_has_filter(h, v, const_cast<char*>("overlap_vntr"))==1) 
	n_filter |= FILTER_MASK_OVERLAP_VNTR;            
    }
    else if (vtype==VT_VNTR) {
      rid = bcf_get_rid(v);
      beg1 = bcf_get_pos1(v) - 1;
      end1 = bcf_get_end_pos1(v) + 1;
      
      char *motif = NULL;
      int32_t n = 0;
      
      if (bcf_get_info_string(h, v, "MOTIF", &motif, &n)>0) {
	this->motif.assign(motif);
	free(motif);
      }
    }
    else {
      return;
    }

    //pls = NULL; //(uint8_t*)calloc( nsamples*3, sizeof(uint8_t) );
    //ads = NULL; //(uint8_t*)calloc( nsamples*3, sizeof(uint8_t) );
}

/**
 * Clears this record.
 */
void JointGenotypeBlockRecord::clearTemp() {
  tmp_dp_q20 = 0;
  tmp_dp_ra = 0;
  tmp_bq_s1 = tmp_bq_s2 = 0;
  tmp_mq_s1 = tmp_mq_s2 = 0;
  tmp_cy_s1 = tmp_cy_s2 = 0;
  tmp_st_s1 = tmp_st_s2 = 0;
  tmp_al_s1 = tmp_bq_al = tmp_mq_al = tmp_cy_al = tmp_st_al = tmp_nm_al = 0;
  tmp_nm_s1 = tmp_nm_s2 = 0;
  tmp_oth_exp_q20 = tmp_oth_obs_q20 = 0;
  tmp_pls[0] = tmp_pls[1] = tmp_pls[2] = 1.;
  tmp_ads[0] = tmp_ads[1] = tmp_ads[2] = 0;
}

void JointGenotypeBlockRecord::clear()
{
  //v = NULL;
  vtype = -1;

  //no_nonref = 0;

  //pls = ads = NULL;
  
  bqr_num = bqr_den = 0;
  mqr_num = mqr_den = 0;
  cyr_num = cyr_den = 0;
  str_num = str_den = 0;
  nmr_num = nmr_den = 0;
  ior_num = ior_den = 0;
  nm0_num = nm0_den = 0;
  nm1_num = nm1_den = 0;
  abe_num = abe_den = 0;
  abz_num = abz_den = 0;
  //ns_nref = dp_sum = max_gq = 0;
  clearTemp();
}

/**
 * Destructor.
 */
JointGenotypeBlockRecord::~JointGenotypeBlockRecord()
{
  //if (v) bcf_destroy(v);
  //if ( pls ) free(pls);
  //if ( ads ) free(ads);
    if ( alleles.l > 0 ) free(alleles.s);
    //if ( est ) delete est;
}

bcf1_t* JointGenotypeBlockRecord::flush_variant(bcf_hdr_t* hdr, sex_ploidy_map& spmap, uint8_t* pls, uint8_t* ads, frequency_estimator* pFreqEst) {
  
  bcf1_t *nv = bcf_init();
  bcf_clear(nv);
  bcf_set_n_sample(nv, nsamples);

  bcf_set_rid(nv, rid);
  bcf_set_pos1(nv, pos1);
  bcf_update_alleles_str(hdr, nv, alleles.s);

  int32_t* gt = (int32_t*) calloc ( nsamples * 2, sizeof(int32_t) );
  int32_t* pl = (int32_t*) calloc ( nsamples * 3, sizeof(int32_t) );
  int32_t* ad = (int32_t*) calloc ( nsamples * 2, sizeof(int32_t) );
  int32_t* td = (int32_t*) calloc ( nsamples * 1, sizeof(int32_t) );
  int32_t* gq = (int32_t*) calloc ( nsamples * 1, sizeof(int32_t) );
  //float MLE_HWE_AF[2];
  //float MLE_HWE_GF[3];
  double gp = 0, gp_sum = 0, max_gp = 0;
  int32_t best_gt = 0;
  int32_t best_a1 = 0, best_a2 = 0;
  int32_t an = 0;
  int32_t acs[2] = {0,0};
  int32_t gcs[3] = {0,0,0};
  float afs[3];
  int32_t max_gq = 0;
  int64_t dp_sum = 0;

  // obtain PL
  std::copy( pls, pls + (nsamples*3), pl );
  for(int32_t i=0; i < nsamples; ++i) {
    dp_sum += ( td[i] = ( (ad[2*i] = ads[3*i]) + (ad[2*i+1] = ads[3*i+1]) + ads[3*i+2] ) );
  }

  //notice("flush_variant pos1=%d, pl=(%d,%d,%d,%d,%d,%d)",pos1,pl[0],pl[1],pl[2],pl[3],pl[4],pl[5]);

  //int32_t n = 0;
  int32_t adSumHet[2] = {0,0};

  // calculate allele frequency
  int8_t* ploidies = spmap.get_ploidies(nv);
  pFreqEst->set_variant(nv, ploidies, pl);
  pFreqEst->estimate_isaf_em(); // compute pooled_allele frequency
  pFreqEst->score_test_hwe(true);

  //bool isX = ( spmap.get_ploidy_type(nv) == PLOIDY_TYPE_X );
  float* ifs = pFreqEst->ifs;

  //fprintf(stderr,"%d %d\n", (int32_t)isX, spmap.vSex.size());

  for(int32_t i=0; i < nsamples; ++i) {
    int32_t* pli = &pl[ i * 3 ];

    if ( ploidies[i] == 1 ) {
      max_gp = gp_sum = gp = ( phredConv.toProb((uint32_t)pli[0]) * (1.0 - ifs[i]) );
      best_gt = 0; best_a1 = 0; best_a2 = 0;
      gp = ( phredConv.toProb((uint32_t)pli[2]) * ifs[i] );
      gp_sum += gp;
      if ( max_gp < gp ) {
	max_gp = gp;
	best_gt = 2; best_a1 = 1; best_a2 = 1;
      }	
    }
    else if ( ploidies[i] == 2 ) {
      max_gp = gp_sum = gp = ( phredConv.toProb((uint32_t)pli[0]) * (1.0-ifs[i]) * (1.0-ifs[i]) );
      best_gt = 0; best_a1 = 0; best_a2 = 0;

      gp = phredConv.toProb((uint32_t)pli[1]) * 2.0 * ifs[i] * (1.0-ifs[i]);
      gp_sum += gp;
      if ( max_gp < gp ) { max_gp = gp; best_gt = 1; best_a1 = 0; best_a2 = 1; }

      gp = phredConv.toProb((uint32_t)pli[2]) * ifs[i] * ifs[i];
      gp_sum += gp;
      if ( max_gp < gp ) { max_gp = gp; best_gt = 2; best_a1 = 1; best_a2 = 1; }      
    }
    else if ( ploidies[i] == 0 ) { // there is no best-guess genotype
      best_gt = 0;
      max_gp = 0;
      gp_sum = 1e-100;
    }
    else 
      error("[E:%s:%d %s] Unexpected ploidy %d", __FILE__, __LINE__, __PRETTY_FUNCTION__, (int32_t)ploidies[i]);
    
    if ( ( ploidies[i] == 2 ) && ( best_gt == 1 ) ) {
      adSumHet[0] += ad[2*i];
      adSumHet[1] += ad[2*i+1];
    }
      
    double prob = 1.-max_gp/gp_sum;  // to calculate GQ
    if ( prob <= 3.162278e-26 )
      prob = 3.162278e-26;
    if ( prob > 1 )
      prob = 1;
    
    gq[i] = (int32_t)phredConv.err2Phred((double)prob);
    
    if ( ( best_gt > 0 ) && ( max_gq < gq[i] ) ) {
      max_gq = gq[i];
    }
    
    gt[2*i]   = ((best_a1 + 1) << 1);
    gt[2*i+1] = ((best_a2 + 1) << 1);	    
    an += 2;             // still use diploid representation of chrX for now.
    ++acs[best_a1];
    ++acs[best_a2];
    ++gcs[best_gt];
  }
    
  for(size_t i=0; i < 2; ++i) {
    afs[i] = acs[i]/(float)an;
  }

  bcf_update_format_int32(hdr, nv, "GT", gt, nsamples * 2);
  bcf_update_format_int32(hdr, nv, "GQ", gq, nsamples );	  
  bcf_update_format_int32(hdr, nv, "AD", ad, nsamples * 2);
  bcf_update_format_int32(hdr, nv, "DP", td, nsamples );	  
  bcf_update_format_int32(hdr, nv, "PL", pl, nsamples * 3);

  float avgdp = (float)dp_sum / (float)nsamples;
  
  nv->qual = (float)max_gq;

  //if ( acs[1] > 0 ) notice("AC=%d, max-gq=%d, QUAL=%f",acs[1], max_gq, nv->qual);
  
  float flt20[21];
  if ( printTmpInfo ) {
    flt20[0] = bqr_num;
    flt20[1] = bqr_den;
    flt20[2] = mqr_num;
    flt20[3] = mqr_den;
    flt20[4] = cyr_num;
    flt20[5] = cyr_den;
    flt20[6] = str_num;
    flt20[7] = str_den;
    flt20[8] = nmr_num;
    flt20[9] = nmr_den;
    flt20[10] = ior_num;
    flt20[11] = ior_den;
    flt20[12] = nm0_num;
    flt20[13] = nm0_den;
    flt20[14] = nm1_num;
    flt20[15] = nm1_den;
    flt20[16] = abe_num;
    flt20[17] = abe_den;
    flt20[18] = abz_num;
    flt20[19] = abz_den;
  }

  bcf_update_info_float(hdr, nv, "AVGDP", &avgdp, 1);	  
  bcf_update_info_int32(hdr, nv, "AC", &acs[1], 1);
  bcf_update_info_int32(hdr, nv, "AN", &an, 1);
  bcf_update_info_float(hdr, nv, "AF", &afs[1], 1);
  bcf_update_info_int32(hdr, nv, "GC", gcs, 3);
  bcf_update_info_int32(hdr, nv, "GN", &nsamples, 1);

  float pooled_af = pFreqEst->pooled_af;
  float hweslp0 = (float)((pFreqEst->hwe0z > 0 ? -1 : 1) * log10( erfc(fabs(pFreqEst->hwe0z)/sqrt(2.0)) + 1e-100 ));
  float hweslp1 = (float)((pFreqEst->hwe1z > 0 ? -1 : 1) * log10( erfc(fabs(pFreqEst->hwe1z)/sqrt(2.0)) + 1e-100 ));  
  
  bcf_update_info_float(hdr, nv, "HWEAF_P", &pooled_af, 1);
  bcf_update_info_float(hdr, nv, "FIBC_P", &(pFreqEst->ibc0), 1);
  bcf_update_info_float(hdr, nv, "HWE_SLP_P", &hweslp0, 1);
  bcf_update_info_float(hdr, nv, "FIBC_I", &(pFreqEst->ibc1), 1);
  bcf_update_info_float(hdr, nv, "HWE_SLP_I", &hweslp1, 1);

  float max_if = 0, min_if = 1;
  for(int32_t j=0; j < nsamples; ++j) {
    if ( ifs[j] > max_if ) max_if = ifs[j];
    if ( ifs[j] < min_if ) min_if = ifs[j];
  }
  bcf_update_info_float(hdr, nv, "MAX_IF", &max_if, 1);
  bcf_update_info_float(hdr, nv, "MIN_IF", &min_if, 1);
  bcf_update_info_float(hdr, nv, "BETA_IF", pFreqEst->betas, pFreqEst->ndims);

  float abe = (adSumHet[0] + 0.5)/(adSumHet[0] + adSumHet[1] + 1.0);
  bcf_update_info_float(hdr, nv, "ABE", &abe, 1);

  // update filter
  if ( n_filter & FILTER_MASK_OVERLAP_SNP )
    bcf_add_filter(hdr, nv, bcf_hdr_id2int(hdr, BCF_DT_ID, "overlap_snp"));
  if ( n_filter & FILTER_MASK_OVERLAP_INDEL )
    bcf_add_filter(hdr, nv, bcf_hdr_id2int(hdr, BCF_DT_ID, "overlap_indel"));
  if ( n_filter & FILTER_MASK_OVERLAP_VNTR )
    bcf_add_filter(hdr, nv, bcf_hdr_id2int(hdr, BCF_DT_ID, "overlap_vntr"));  

  //bcf_update_info_int32(odw->hdr, nv, "NS_NREF", &v_ns_nrefs[k], 1);
  abe_num /= (abe_den+1e-6); bcf_update_info_float(hdr, nv, "ABE",  &abe_num, 1);
  abz_num /= sqrt(abz_den+1e-6); bcf_update_info_float(hdr, nv, "ABZ",  &abz_num, 1);	  	  
  bqr_num /= sqrt(bqr_den+1e-6); bcf_update_info_float(hdr, nv, "BQZ", &bqr_num, 1);	  
  mqr_num /= sqrt(mqr_den+1e-6); bcf_update_info_float(hdr, nv, "MQZ", &mqr_num, 1);	  
  cyr_num /= sqrt(cyr_den+1e-6); bcf_update_info_float(hdr, nv, "CYZ", &cyr_num, 1);	  
  str_num /= sqrt(str_den+1e-6); bcf_update_info_float(hdr, nv, "STZ", &str_num, 1);
  nmr_num /= sqrt(nmr_den+1e-6); bcf_update_info_float(hdr, nv, "NMZ", &nmr_num, 1);
  ior_num = log(ior_num/ior_den+1e-6)/log(10.); bcf_update_info_float(hdr, nv, "IOR", &ior_num, 1);
  nm1_num /= (nm1_den+1e-6); bcf_update_info_float(hdr, nv, "NM1", &nm1_num, 1);
  nm0_num /= (nm0_den+1e-6); bcf_update_info_float(hdr, nv, "NM0", &nm0_num, 1);

  if ( printTmpInfo ) {
    bcf_update_info_float(hdr, nv, "FLT20", flt20, 20);
  }

  //odw->write(nv);
  //bcf_destroy(nv);
  free(gt);
  free(gq);
  free(pl);
  free(ad);
  free(td);
  
  //free(pls); pls = NULL;
  //free(ads); ads = NULL;

  //delete est;
  if ( alleles.l > 0 ) {
    alleles.l = 0;
    free(alleles.s);
  }

  return nv;
}

void JointGenotypeBlockRecord::flush_sample( int32_t sampleIndex, gzFile fp ) {
  //uint8_t* p_pls = &pls[sampleIndex*3];
  //uint8_t* p_ads = &ads[sampleIndex*3];
  uint8_t p_pl_ads[6] = {0,0,0,0,0,0};
  uint8_t* p_pls = &p_pl_ads[0];
  uint8_t* p_ads = &p_pl_ads[3];

  for(int32_t i=0; i < 3; ++i)
    if ( tmp_pls[i] < 1e-300 ) tmp_pls[i] = 1e-300;

  int32_t imax = ( tmp_pls[0] > tmp_pls[1] ) ? ( tmp_pls[0] > tmp_pls[2] ? 0 : 2 ) : ( tmp_pls[1] > tmp_pls[2] ? 1 : 2);
  for(int32_t i=0; i < 3; ++i) {
    uint32_t l = phredConv.err2Phred((double)(tmp_pls[i]/tmp_pls[imax]));
    p_pls[i] = ((l > 255) ? 255 : (uint8_t)l);
    p_ads[i] = ((tmp_ads[i] > 255) ? 255 : (uint8_t)tmp_ads[i]);
  }

  //error("tmp_pls = (%lg, %lg, %lg), p_pls = (%u, %u, %u)", tmp_pls[0], tmp_pls[1], tmp_pls[2], p_pls[0], p_pls[1], p_pls[2]); 

  if ( ( p_pls[0] == 0 ) && ( p_pls[2] == 0 ) && ( p_pls[1] > 0 ) ) 
    error("%le %le %le", tmp_pls[0], tmp_pls[1], tmp_pls[2]);

  float sqrt_dp_ra = sqrt((float)tmp_dp_ra);
  float ior = (float)(tmp_oth_obs_q20 / (tmp_oth_exp_q20 + 1e-6));
  float nm1 = tmp_al_s1 == 0 ? 0 : tmp_nm_al / (float)tmp_al_s1;
  float nm0 = (tmp_dp_ra - tmp_al_s1) == 0 ? 0 : (tmp_nm_s1-tmp_nm_al) / (float)(tmp_dp_ra - tmp_al_s1);
  float w_dp_ra  = log(tmp_dp_ra+1.); //sqrt(dp_ra);
  float w_dp_q20 = log(tmp_dp_q20+1.); //sqrt(dp_q20);
  float w_al_s1  = log(tmp_al_s1+1.); //sqrt(al_s1);
  float w_ref_s1 = log(tmp_dp_ra - tmp_al_s1+1.);
  
  if ( p_pls[1] == 0 ) { // het genotypes
    abe_num += (w_dp_ra * (tmp_dp_ra - tmp_al_s1 + 0.05) / (double)(tmp_dp_ra + 0.1));
    abe_den += w_dp_ra;
    
    // E(r) = 0.5(r+a) V(r) = 0.25(r+a)
    abz_num += w_dp_ra * (tmp_dp_ra - tmp_al_s1 - tmp_dp_ra*0.5)/sqrt(0.25 * tmp_dp_ra + 1e-3);
    abz_den += (w_dp_ra * w_dp_ra);
    
    float bqr = sqrt_dp_ra * compute_correlation( tmp_dp_ra, tmp_bq_al, tmp_bq_s1, tmp_bq_s2, tmp_al_s1, tmp_al_s1, .1 );
    float mqr = sqrt_dp_ra * compute_correlation( tmp_dp_ra, tmp_mq_al, tmp_mq_s1, tmp_mq_s2, tmp_al_s1, tmp_al_s1, .1 );
    float cyr = sqrt_dp_ra * compute_correlation_f( tmp_dp_ra, tmp_cy_al, tmp_cy_s1, tmp_cy_s2, (float)tmp_al_s1, (float)tmp_al_s1, .1 );
    float str = sqrt_dp_ra * compute_correlation( tmp_dp_ra, tmp_st_al, tmp_st_s1, tmp_st_s1, tmp_al_s1, tmp_al_s1, .1 );
    float nmr = sqrt_dp_ra * compute_correlation( tmp_dp_ra, tmp_nm_al, tmp_nm_s1, tmp_nm_s2, tmp_al_s1, tmp_al_s1, .1 );
    
    // Use Stouffer's method to combine the z-scores, but weighted by log of sample size
    bqr_num += (bqr * w_dp_ra); bqr_den += (w_dp_ra * w_dp_ra);
    mqr_num += (mqr * w_dp_ra); mqr_den += (w_dp_ra * w_dp_ra);
    cyr_num += (cyr * w_dp_ra); cyr_den += (w_dp_ra * w_dp_ra);
    str_num += (str * w_dp_ra); str_den += (w_dp_ra * w_dp_ra);
    nmr_num += (nmr * w_dp_ra); nmr_den += (w_dp_ra * w_dp_ra);	  
  }
  
  ior_num += (ior * w_dp_q20); ior_den += w_dp_q20;
  nm1_num += (nm1 * w_al_s1);  nm1_den += w_al_s1;
  nm0_num += (nm0 * w_ref_s1); nm0_den += w_ref_s1;

  //notice("sampleIndex = %d, pos1 = %d, fp = %x", sampleIndex, pos1, fp);
  //notice("Flushing sample %d, pos1 %d, fp %x, fp->fn %s, fp->fp.bgzf %x, fp->fp.hfile %x, %d %d %d %d %d", sampleIndex, pos1, fp, fp->fn, fp->fp.bgzf, fp->fp.hfile, fp->is_bin, fp->is_write, fp->is_be, fp->is_cram, fp->is_bgzf);                                   


  //if ( bgzf_write(fp->fp.bgzf, p_pl_ads, sizeof(uint8_t)*6) < 0 )
  if ( gzwrite(fp, p_pl_ads, sizeof(uint8_t)*6) != 6 )
    error("Error in writing sample %d at pos %d", sampleIndex, pos1);

  clearTemp();
}

void JointGenotypeBlockRecord::add_allele( double contam, int32_t allele, uint8_t mapq, bool fwd, uint32_t q, int32_t cycle, uint32_t nm ) {
  //if ( pos1 == 10000004 ) notice("add_allele(%lg,%d,%u,%d,%u,%d,%u) called",contam,allele,mapq,fwd ? 0 : 1,q,cycle,nm);
  if ( q > 40 )
    q = 40;

  double pe = phredConv.toProb((uint32_t)q);
  if ( pe > 0.75 ) pe = 0.75;
  double pm = 1 - pe;
  
  if ( allele == 0 ) {
    ++tmp_ads[0];
    tmp_pls[0] *= ( pm * (1.-contam) + pe * contam / 3. );
    tmp_pls[1] *= ( pm / 2. + pe / 6. );
    tmp_pls[2] *= ( pm * contam + pe * (1.-contam) / 3. );
  }
  else if ( allele > 0 ) { // currently, bi-allelic only
    ++tmp_ads[1];
    tmp_pls[0] *= ( pm * contam + pe * (1.-contam) / 3. );
    tmp_pls[1] *= ( pm / 2. + pe / 6. );
    tmp_pls[2] *= ( pm * (1.-contam) + pe * contam / 3. );    
  }
  else {
    ++tmp_ads[2];
  }
  double sump = tmp_pls[0] + tmp_pls[1] + tmp_pls[2];
  if ( sump < 1e-300 ) sump = 1e-300;
  tmp_pls[0] /= sump;
  tmp_pls[1] /= sump;
  tmp_pls[2] /= sump;

  if ( allele >= 0 ) {
    if ( q > 20 ) {
      tmp_oth_exp_q20 += (phredConv.toProb((uint32_t)q) * 2. / 3.);
      ++tmp_dp_q20;
    }

    float log_td = (cycle > 0) ? 0-logf((float)cycle) : 0;

    ++tmp_dp_ra;
    tmp_bq_s1 += q;
    tmp_bq_s2 += (q*q);
    tmp_mq_s1 += mapq;
    tmp_mq_s2 += (mapq * mapq);
    tmp_cy_s1 += log_td;
    tmp_cy_s2 += (log_td * log_td);
    tmp_st_s1 += fwd;

    if ( allele > 0 ) {
      ++tmp_al_s1;
      tmp_bq_al += q;
      tmp_mq_al += mapq;
      tmp_cy_al += log_td;
      tmp_st_al += fwd;
      tmp_nm_al += (nm-1);
      tmp_nm_s1 += (nm-1);
      tmp_nm_s2 += (nm-1)*(nm-1);
    }
    else {
      tmp_nm_s1 += nm;
      tmp_nm_s2 += (nm * nm);
    }
  }
  else {
    if ( q > 20 ) {
      tmp_oth_exp_q20 += (phredConv.toProb(q) * 2. / 3.);
      ++tmp_oth_obs_q20;
      ++tmp_dp_q20;
    }
  }
}

/**
 * Collects sufficient statistics from read for variants to be genotyped.
 */
void JointGenotypeBlockRecord::process_read(AugmentedBAMRecord& as, int32_t sampleIndex, double contam) {
  if (vtype==VT_SNP) {
    if (v_alleles.size()==2) {
      bam1_t *s = as.s;
      
      char strand = bam_is_rev(s) ? 'R' : 'F';  
      int32_t allele = 0;
      //uint32_t bpos1 = bam_get_pos1(s);
      //uint8_t* seq = bam_get_seq(s);
      uint8_t* qual = bam_get_qual(s);
      int32_t rlen = bam_get_l_qseq(s);
      uint8_t mapq = bam_get_mapq(s);
      uint32_t q = 30;
      int32_t cycle = 0;

      std::vector<uint32_t>& aug_cigar = as.aug_cigar;
      //std::vector<std::string>& aug_ref = as.aug_ref;
      std::vector<std::string>& aug_alt = as.aug_alt;
      
      int32_t vpos1 = pos1;
      int32_t cpos1 = bam_get_pos1(s);
      int32_t rpos0 = 0;

      for (uint32_t i=0; i<aug_cigar.size(); ++i) {
	uint32_t oplen = bam_cigar_oplen(aug_cigar[i]);
	char opchr = bam_cigar_opchr(aug_cigar[i]);

	if (opchr=='S') {
	  rpos0 += oplen;
	}
	else if ( opchr=='H') {} 
	else if (opchr=='=') {
	  if (vpos1>=cpos1 && vpos1 <= (int32_t)(cpos1+oplen-1)) {
	    rpos0 += vpos1-cpos1;
	    allele = 0;
	    q = qual[rpos0];
	    cycle = rpos0<(rlen>>1) ? (rpos0+1) : -(rlen - rpos0 + 1);
	    break;
	  }
	  else {
	    cpos1 += oplen;
	    rpos0 += oplen;
	  }
	}
	else if (opchr=='X') {
	  if (vpos1==cpos1) {
	    allele = (aug_alt[i].at(0) == v_alleles[1].at(0)) ? 1 : -1;
	    q = qual[rpos0];
	    cycle = rpos0<(rlen>>1) ? (rpos0+1) : -(rlen - rpos0 + 1);
	    break;
	  }
	  
	  ++cpos1;
	  ++rpos0;
	}
	else if (opchr=='I') {
	  rpos0 += oplen;
	}
	else if (opchr=='D') {
	  cpos1 += oplen;
	}
	else {
	  std::cerr << "unrecognized cigar state " << opchr << "\n";
	  //            exit(1);
	}
      }

      uint32_t no_mismatches = as.no_mismatches;
      if (allele!=0 && q<20) {
	++no_mismatches;
      }
      if (allele!=0 && no_mismatches==0) {
	//no_mismatches = 1;
	std::cerr << "something wrong\n";
      }
      
      add_allele( contam, allele, mapq, strand == 'F', q, cycle, no_mismatches );
    }
    else { //multiallelic
      //abort();
    }
  }
  else if (vtype==VT_INDEL) {
    if (v_alleles.size()==2) {
      if (as.beg1 <= beg1 && end1 <= as.end1) {
	bam1_t *s = as.s;
	
	char strand = bam_is_rev(s) ? 'R' : 'F';
	int32_t allele = 0;
	//uint32_t bpos1 = bam_get_pos1(s);
	//uint8_t* seq = bam_get_seq(s);
	//uint8_t* qual = bam_get_qual(s);
	uint32_t rlen = bam_get_l_qseq(s);
	uint8_t mapq = bam_get_mapq(s);
	
	uint32_t q = len*30;
	uint32_t cycle = 10;
	
	std::vector<uint32_t>& aug_cigar = as.aug_cigar;
	std::vector<std::string>& aug_ref = as.aug_ref;
	std::vector<std::string>& aug_alt = as.aug_alt;
	
	int32_t vpos1 = pos1;
	
	int32_t cpos1 = bam_get_pos1(s);
	int32_t rpos0 = 0;
	
	for (uint32_t i=0; i<aug_cigar.size(); ++i) {
	  uint32_t oplen = bam_cigar_oplen(aug_cigar[i]);
	  char opchr = bam_cigar_opchr(aug_cigar[i]);
	  
	  if (opchr=='S') {
	    rpos0 += oplen;
	  }
	  else if (opchr=='H') {}
	  else if (opchr=='=') {
	    if (vpos1>=cpos1 && vpos1<=(int32_t)(cpos1+oplen-1)) {
	      cycle = strand == 'F' ? (rpos0+1) : (rlen - rpos0);
	    }
	    cpos1 += oplen;
	    rpos0 += oplen;
	  }
	  else if (opchr=='X') {
	    if (cpos1-1==vpos1) {
	      cycle = strand == 'F' ? (rpos0+1) : (rlen - rpos0);
	    }
	    ++cpos1;
	    ++rpos0;
	  }
	  else if (opchr=='I') {
	    if (dlen>0 && cpos1-1==vpos1) {
	      if (indel==aug_alt[i]) {
		q = len*30;
		allele = 1;
	      }
	      else {
		q = abs(len-aug_ref[i].size())*30;
		allele = -1;
	      }
	      
	      cycle = strand == 'F' ? (rpos0+1) : (rlen - rpos0);
	      break;
	    }
	    else if (dlen<0 && cpos1-1==vpos1) {
	      q = 30;
	      allele = -3;
	      cycle = strand == 'F' ? (rpos0+1) : (rlen - rpos0);
	      break;
	    }

	    rpos0 += oplen;
	  }
	  else if (opchr=='D') {
	    if (dlen<0 && cpos1-1==vpos1) {
	      if (indel==aug_ref[i]) {
		q = len*30;
		allele = 1;
	      }
	      else {
		q = abs(len-aug_ref[i].size())*30;
		allele = -1;
	      }
	      cycle = strand == 'F' ? (rpos0+1) : (rlen - rpos0);
	      break;
	    }
	    else if (dlen>0 && cpos1-1==vpos1) {
	      q = 30;
	      allele = -2;
	      
	      cycle = strand == 'F' ? (rpos0+1) : (rlen - rpos0);
	      break;
	    }
	    
	  }
	  else {
	    std::cerr << "unrecognized cigar state " << opchr << "\n";
	  }
	}

	add_allele( contam, allele, mapq, strand == 'F', q, rand() % 75, as.no_mismatches );		
      }
      else {
	bam1_t *s = as.s;
	uint8_t mapq = bam_get_mapq(s);

	add_allele( contam, -1, mapq, bam_is_rev(as.s) ? false : true, 20, rand() % 75, as.no_mismatches );		
      }
    }
    else { //multiallelic
      //abort();
    }
  }
  else if (vtype==VT_VNTR) {
    //abort();
 //    bam1_t *s = as.s;
    
//     char strand = bam_is_rev(s) ? 'R' : 'F';
//     int32_t allele = 0;
//     uint32_t pos1 = bam_get_pos1(s);
//     uint8_t* seq = bam_get_seq(s);
//     uint8_t* qual = bam_get_qual(s);
//     uint32_t rlen = bam_get_l_qseq(s);
//     uint8_t mapq = bam_get_mapq(s);
    
//     std::vector<uint32_t>& aug_cigar = as.aug_cigar;
//     std::vector<std::string>& aug_ref = as.aug_ref;
//     std::vector<std::string>& aug_alt = as.aug_alt;
    
//     //genomic bookend positions of VNTR
//     int32_t vpos1 = g->beg1-1;
//     int32_t vend1 = g->end1+1;

//     //position with respect to read
//     int32_t cpos1 = bam_get_pos1(s);
//     int32_t rpos0 = 0;
    
//     //genomic bookend positions of VNTR translated to read position
//     int32_t pos0 = -1;
//     int32_t end0 = -1;
    
    
//     //locate repeat region on read.
//     //translating genomic coordinates to read positions
//     for (uint32_t i=0; i<aug_cigar.size(); ++i) {
//       uint32_t oplen = bam_cigar_oplen(aug_cigar[i]);
//       char opchr = bam_cigar_opchr(aug_cigar[i]);
      
//       if (opchr=='S') {
// 	rpos0 += oplen;
//       }
//       else if (opchr=='=') {
// 	if (pos0==-1 && (cpos1<=vpos1 && vpos1<=(cpos1+oplen-1))) {
// 	  pos0 = rpos0 + (vpos1-cpos1+1);
// 	}
	
// 	if (end0==-1 && (cpos1<=vend1 && vend1<=(cpos1+oplen-1))) {
// 	  end0 = rpos0 + (vend1-cpos1+1);
// 	  break;
// 	}

// 	cpos1 += oplen;
// 	rpos0 += oplen;
//       }
//       else if (opchr=='X') {
// 	if (pos0==-1 && (cpos1==vpos1)) {
// 	  pos0 = rpos0;
// 	}
	
// 	if (end0==-1 && (cpos1==vend1)) {
// 	  end0 = rpos0;
// 	  break;
// 	}
	
// 	++cpos1;
// 	++rpos0;
//       }
//       else if (opchr=='I') {
// 	rpos0 += oplen;
//       }
//       else if (opchr=='D') {
// 	cpos1 += oplen;
//       }
//       else {
// 	std::cerr << "unrecognized cigar state " << opchr << "\n";
//       }
//     }
    
//     //compute repeat tract units
//     float counts = 0;
    
// //        std::cerr << "pos0,end0 = " << pos0 << "," <<end0 <<  " (" << g->motif.size() <<")\n";

//     if (pos0!=-1 && end0!=-1) {
//       counts = ((float)(end0-pos0+1))/g->motif.size();
//     }
    
//     if (counts) {
//       //update genotyping record
//       ++g->depth;
//       g->counts.push_back(counts);
//       g->strands.append(1, strand);
//       g->no_mismatches.push_back(as.no_mismatches);
  }
}
