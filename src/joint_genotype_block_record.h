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

#ifndef JOINT_GENOTYPE_BLOCK_RECORD_H
#define JOINT_GENOTYPE_BLOCK_RECORD_H

#include <zlib.h>
#include "htslib/vcf.h"
#include "htslib/faidx.h"
#include "bcf_ordered_writer.h"
#include "variant.h"
#include "hts_utils.h"
#include "augmented_bam_record.h"
#include "sex_ploidy_map.h"
#include "frequency_estimator.h"
#include "Eigen/Dense"

#define FILTER_MASK_OVERLAP_SNP   0x0001
#define FILTER_MASK_OVERLAP_INDEL 0x0002
#define FILTER_MASK_OVERLAP_VNTR  0x0004

/**
 * A generic record that holds information for genotyping a
 * variant across multiple samples.
 *
 * Maintains read information and allows for additional reads
 * till VCF record can be printed out.
 */
class JointGenotypeBlockRecord
{
    public:
    bcf_hdr_t *h;
    //bcf1_t *v;
    int32_t rid;
    int32_t pos1; //position of variant
    //[beg1,end1] is the required overlapping of the variant against the aligned read necessary to make a genotype call.
    //for SNPs, beg1=end1=pos1
    //
    //for Indels, this refers to the flanking positions
    //   insertion 
    //        if T/TG - beg1=pos1, end1=pos1+1
    //        if T/GT - beg1=pos1-1, end1=pos1
    //   deletion  
    //        if TG/T - beg1=pos1, end1=pos1+length(REF)
    //        if TG/G - beg1=pos1-1, end1=pos1+length(REF)-1
    int32_t beg1;
    int32_t end1;
    int32_t vtype;

    //indel specific record
    int32_t dlen;
    uint32_t len;
    std::string indel;

    //vntr specific record
    std::string motif;

    //vntr specific record
    //std::vector<float> counts;

    // sample level information
    int32_t nsamples;
    kstring_t alleles;
    std::vector<std::string> v_alleles;
    uint32_t n_filter;

    //uint8_t* pls;
    //uint8_t* ads;
    //Estimator* est;

    // sufficient statistics for computing INFO field
    float bqr_num, bqr_den; // TMP_0, TMP_1
    float mqr_num, mqr_den; // TMP_2, TMP_3
    float cyr_num, cyr_den; // TMP_4, TMP_5
    float str_num, str_den; // TMP_6, TMP_7
    float nmr_num, nmr_den; // TMP_8, TMP_9
    float ior_num, ior_den; // TMP_10, TMP_11
    float nm0_num, nm0_den; // TMP_12, TMP_13
    float nm1_num, nm1_den; // TMP_14, TMP_15
    float abe_num, abe_den; // TMP_16, TMP_17
    float abz_num, abz_den; // TMP_18, TMP_19
    float tmp_cy_s1, tmp_cy_s2; // TMP_20, TMP_21
    float tmp_cy_al;       // TMP_22
    float tmp_oth_exp_q20, tmp_oth_obs_q20; // TMP_23, TMP_24    
    //float ns_nref, dp_sum, max_gq;

    int32_t tmp_dp_q20;           // TMP_1
    int32_t tmp_dp_ra;            // TMP_2
    int32_t tmp_bq_s1, tmp_bq_s2; // TMP_3, TMP_4
    int32_t tmp_mq_s1, tmp_mq_s2; // TMP_5, TMP_6
    int32_t tmp_st_s1, tmp_st_s2; // TMP_7, TMP_8
    int32_t tmp_al_s1, tmp_bq_al, tmp_mq_al; // TMP_9, TMP_10
    int32_t tmp_st_al, tmp_nm_al; // TMP_11, TMP_12
    int32_t tmp_nm_s1, tmp_nm_s2; // TMP_13, TMP_14
    
    double tmp_pls[3];
    double tmp_ads[3];

    // temporary information to be cleared out per-sample basis
    bool printTmpInfo;

    /**
     * Constructor.
     * @v - VCF record.
     */
    JointGenotypeBlockRecord(bcf_hdr_t*h, bcf1_t *v, int32_t vtype, int32_t nsamples, bool _printTmpInfo);

    /**
     * Clears this record.
     */
    void clear();
    void clearTemp();
    bcf1_t* flush_variant(bcf_hdr_t* hdr, sex_ploidy_map& spmap, uint8_t* pls, uint8_t* ads, frequency_estimator* pFreqEst);
    //void flush_sample( int32_t sampleIndex, htsFile* fp );
    void flush_sample( int32_t sampleIndex, gzFile fp );
    void add_allele( double contam, int32_t allele, uint8_t mapq, bool fwd, uint32_t q, int32_t cycle, uint32_t nm );
    void process_read(AugmentedBAMRecord& as, int32_t sampleIndex, double contam);    

    /**
     * Destructor.
     */
    ~JointGenotypeBlockRecord();

    float compute_correlation(int32_t n, int32_t xy, int32_t x1, int32_t x2, int32_t y1, int32_t y2, float buffer) {
      if ( n == 0 ) return 0;
      float xsd = x2/(float)n - (x1/(float)n)*(x1/(float)n);
      float ysd = y2/(float)n - (y1/(float)n)*(y1/(float)n);
      return ( ( xy/(float)n - x1 * y1 / (float) n / (float) n ) / sqrt( xsd * ysd + buffer ) );
    }

    float compute_correlation_f(int32_t n, float xy, float x1, float x2, float y1, float y2, float buffer) {
      if ( n == 0 ) return 0;
      float xsd = x2/(float)n - (x1/(float)n)*(x1/(float)n);
      float ysd = y2/(float)n - (y1/(float)n)*(y1/(float)n);
      return ( ( xy/(float)n - x1 * y1 / (float) n / (float) n ) / sqrt( xsd * ysd + buffer ) );
    }
    
};

#endif
