/* The MIT License

   Copyright (c) 2015 Adrian Tan <atks@umich.edu> and Hyun Min Kang <hmkang@umich.edu>

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

#ifndef JOINT_GENOTYPE_BLOCK_READER_H
#define JOINT_GENOTYPE_BLOCK_READER_H

#include <string>
#include <vector>
#include <zlib.h>
#include "htslib/kseq.h"
#include "htslib/vcf.h"
#include "hts_utils.h"
#include "joint_genotype_block_record.h"
#include "bcf_ordered_reader.h"
#include "variant.h"
#include "variant_manip.h"
#include "augmented_bam_record.h"
#include "sex_ploidy_map.h"
#include "frequency_estimator.h"
#include "Eigen/Dense"

/**
 * Wrapper for BCFOrderedReader.
 *
 * VCF records are wrapped in GenotyingRecord and are
 * maintained in a buffer.
 *
 */
class JointGenotypeBlockReader
{
    public:

    ///////
    //i/o//
    ///////
    BCFOrderedReader* odr; // anchor VCF
    int32_t unit;          // Number of variants per block

    //////////////////
    //buffer related//
    //////////////////
    std::vector<JointGenotypeBlockRecord*> gRecords;
    std::string chrom;
    AugmentedBAMRecord as;
    Variant variant;

    int32_t lastFirst;
    //int32_t currentSampleIndex;

    ///////////
    //options//
    ///////////
    bool output_annotations;
    

    /////////
    //stats//
    /////////
    uint32_t no_snps_genotyped;
    uint32_t no_indels_genotyped;
    uint32_t no_vntrs_genotyped;

    /////////
    //tools//
    /////////
    VariantManip *vm;
    //LogTool lt;

    std::string tmp_prefix;
    std::vector<std::string> sample_names;
    std::vector<double> sample_contams;
    Eigen::MatrixXd eV;
    //Eigen::BDCSVD<Eigen::MatrixXd>* pSVD;
    frequency_estimator* pFreqEst;

    std::vector<int32_t> blockStarts;
    std::vector<int32_t> blockEnds;
    std::vector<int32_t>  blockIDs;
    //std::vector<htsFile*> blockFHs;
    std::vector<gzFile> blockFHs;
    std::vector<std::string> blockFNs;
    uint8_t* blockPLs;
    uint8_t* blockADs;
    

    /**
     * Constructor.
     */
    JointGenotypeBlockReader(std::string in_vcf_filename, std::vector<GenomeInterval>& intervals, std::string out_tmp_prefix, int32_t nsamples, int32_t nUnit, bool printTmpInfo);
    ~JointGenotypeBlockReader();

    void set_sample(int32_t sampleIndex, const char* sampleName, double contam, std::vector<double>& evec);
    void flush_sample(int32_t sampleIndex);
    void close_blocks();
    void remove_temp_files();
    
    /**
     * Collects sufficient statistics from read for variants to be genotyped.
     */
    int32_t process_read(bam_hdr_t *h, bam1_t *s, int32_t sampleIndex);

    inline int32_t numVariants() { return (int32_t)gRecords.size(); }
    bcf1_t* flush_variant(int32_t variantIndex, bcf_hdr_t* hdr, sex_ploidy_map& spmap);
    void write_header(BCFOrderedWriter* odw, bool printTmpInfo = false);
};

#endif
