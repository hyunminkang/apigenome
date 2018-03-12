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

#include "joint_genotype_block_reader.h"

/**
 * Constructor.
 */
JointGenotypeBlockReader::JointGenotypeBlockReader(std::string filename, std::vector<GenomeInterval>& intervals, std::string out_tmp_prefix, int32_t nsamples, int32_t nUnit, bool printTmpInfo)
{
    vm = new VariantManip();  

    // read input BCF files and create genotyping records for every variants
    odr = new BCFOrderedReader(filename, intervals);

    tmp_prefix = out_tmp_prefix;
    unit = nUnit;
    
    bcf1_t* v = bcf_init();
    int32_t vidx = 0;
    int32_t bidx = 0;
    blockStarts.push_back(0);
    blockEnds.push_back(-1);
    while( odr->read(v) ) {
      int32_t vtype = vm->classify_variant(odr->hdr, v, variant);
      //if ( (rand() % 10000) == 0 )
      //notice("foo");
      int32_t pos1 = bcf_get_pos1(v);
      
      if ( ( pos1 < intervals[0].start1 ) || ( pos1 > intervals[0].end1 ) ) continue;
	
      if ( ( vtype == VT_VNTR ) || ( bcf_get_n_allele(v) != 2 ) ) continue;

      JointGenotypeBlockRecord* jgr = new JointGenotypeBlockRecord(odr->hdr, v, vtype, nsamples, printTmpInfo);
      gRecords.push_back(jgr);
      blockEnds.back() = vidx;
      blockIDs.push_back(bidx);
      ++vidx;
      if ( vidx % unit == 0 ) {
	++bidx;
	blockStarts.push_back(vidx);
	blockEnds.push_back(vidx-1);
      }
    }
    bcf_destroy(v);

    //odr->close();
    //delete odr;

    //////////////////////////
    //options initialization//
    //////////////////////////
    output_annotations = false;

    ////////////////////////
    //stats initialization//
    ////////////////////////
    no_snps_genotyped = 0;
    no_indels_genotyped = 0;
    no_vntrs_genotyped = 0;

    lastFirst = 0;
    //currentSampleIndex = -1;

    ////////////////////////
    //tools initialization//
    ////////////////////////
    sample_names.resize(nsamples);
    sample_contams.resize(nsamples);

    /// Create temporary files to store binary data - 6 bytes per sample/variants

    if ( blockStarts.back() > blockEnds.back() ) {
      //error("foo");
      --bidx;
      if ( bidx >= 0 ) {
        blockStarts.resize(bidx);
        blockEnds.resize(bidx);
      }
    }
    //error("%d", bidx);

    char buf[255];

    //htsFormat fmt;
    //fmt.category = unknown_category;
    //fmt.format = binary_format;
    //fmt.compression = no_compression;

    for(int32_t i=0; i <= bidx; ++i) {
      sprintf(buf,".%08d.gz",i);
      //sprintf(buf,".%08d.bin",i);
      std::string tmp_file_name = tmp_prefix + buf;
      blockFNs.push_back(tmp_file_name);
      //htsFile* fp = hts_open(tmp_file_name.c_str(), "wb");
      //htsFile* fp = hts_open_format(tmp_file_name.c_str(), "w", &fmt);
      gzFile fp = gzopen(tmp_file_name.c_str(),"wb");
      if ( fp == NULL )
	error("Cannot create temporary file %s", tmp_file_name.c_str());
      blockFHs.push_back(fp);
    }

    notice("# Variants = %d",vidx);

    blockPLs = NULL;
    blockADs = NULL;
    //pSVD = NULL;
    pFreqEst = NULL;
}

JointGenotypeBlockReader::~JointGenotypeBlockReader() {
  //notice("foo");
  delete vm;
  delete odr;
  for(int32_t i=0; i < (int32_t)gRecords.size(); ++i) {
    delete gRecords[i];
  }
  remove_temp_files();
}

void JointGenotypeBlockReader::flush_sample(int32_t sampleIndex) {
  /*
  if ( sampleIndex < 0 )
    sampleIndex = currentSampleIndex;
  */

  for(int32_t i=0; i < (int)gRecords.size(); ++i) {
    //htsFile* fp = blockFHs[blockIDs[i]];
    //gzFile fp = blockFHs[blockIDs[i]];
    //notice("Flushing sample %d, block %d, variant %d / %u", sampleIndex, blockIDs[i], i, gRecords.size());
    //notice("Flushing sample %d, block %d, variant %d, fp %x, fp->fn %s, fp->fp.bgzf %x, fp->fp.hfile %x, %d %d %d %d %d", sampleIndex, blockIDs[i], i, blockFHs[blockIDs[i]], fp->fn, fp->fp.bgzf, fp->fp.hfile, fp->is_bin, fp->is_write, fp->is_be, fp->is_cram, fp->is_bgzf);
    gRecords[i]->flush_sample( sampleIndex, blockFHs[blockIDs[i]] );
  }
}

void JointGenotypeBlockReader::remove_temp_files() {
  if ( blockPLs ) { free(blockPLs); blockPLs = NULL; }
  if ( blockADs ) { free(blockADs); blockADs = NULL; }

  for(int32_t i=0; i < (int32_t)blockFNs.size(); ++i) {
    if ( remove(blockFNs[i].c_str()) ) {
      error("Error in removing temporary file %s", blockFNs[i].c_str());
    }
  }
}


void JointGenotypeBlockReader::close_blocks() {
  for(int32_t i=0; i < (int32_t)blockFHs.size(); ++i) {
    if ( blockFHs[i] == NULL ) continue;
    notice("Closing file %s",blockFNs[i].c_str());
    //if ( hts_close(blockFHs[i]) ) {
    if ( gzclose(blockFHs[i]) ) {
      error("Error in closing block %d", i);
    }
    blockFHs[i] = NULL;
  }
}

void JointGenotypeBlockReader::set_sample(int32_t sampleIndex, const char* sampleName, double contam, std::vector<double>& evec) {
  lastFirst = 0;
  sample_names[sampleIndex] = sampleName;
  sample_contams[sampleIndex] = contam;
  if ( eV.cols() == 0 ) {
    eV = Eigen::MatrixXd::Zero((int32_t)sample_names.size(), evec.size() + 1);
  }
  else if ( eV.cols() != (int32_t) evec.size() + 1 ) {
    error("[E:%s:%d %s] Assertion failed: ev.cols() = %d != 1 + evec.size() = %d", __FILE__, __LINE__, __PRETTY_FUNCTION__, eV.cols(), (int32_t)evec.size());
  }
  eV(sampleIndex,0) = 1;
  for(int32_t i=0; i < (int32_t) evec.size(); ++i ) {
    eV(sampleIndex, i+1) = evec[i];
  }
  //currentSampleIndex = sampleIndex;
}

/**
 * Collects sufficient statistics from read for variants to be genotyped.
 *
 * The VCF records in the buffer must never occur before
 */
int32_t JointGenotypeBlockReader::process_read(bam_hdr_t *h, bam1_t *s, int32_t sampleIndex)
{
    //wrap bam1_t in AugmentBAMRecord
    as.initialize(h, s);

    int32_t tid = (int32_t)bam_get_tid(s);
    int32_t beg1 = (int32_t)as.beg1;
    int32_t end1 = (int32_t)as.end1;

    int32_t nvisited = 0;

    //collect statistics for variant records that are in the buffer and overlap with the read
    JointGenotypeBlockRecord* jgr;
    for(int32_t i = lastFirst; i < (int)gRecords.size(); ++i) {
      jgr = gRecords[i];
      
      //same chromosome
      if (tid == jgr->rid) {
	if (end1 < jgr->beg1) // read ends earlier than the last record to visit -- no need to investigate
	  return nvisited;
	else if (beg1 > jgr->end1) { // read begins later than the last record to visit -- advance the last record
	  ++lastFirst;
	  continue;
	}
	else if ((beg1 <= jgr->pos1) && (jgr->pos1 <= end1)) { // variant position is covered by the read
	  ++nvisited;
	  //jgr->process_read(as, currentSampleIndex, sample_contams[currentSampleIndex]);
	  jgr->process_read(as, sampleIndex, sample_contams[sampleIndex]);
	  //	  if (beg1 <= jgr->beg1 && jgr->end1 <= end1) {
	  //                    std::cerr << "COMPLETE";
	  //                }
	  //                else
	  //                {
	  //drop
	  //                    std::cerr << "PARTIAL";
          //      }
          //  }
          //  else
          //  {
          //      //other types of overlap, just ignore
          //  }
	  //            std::cerr << "\n";
        }
      }
      else if ( tid < jgr->rid )
	return nvisited;
      else if ( tid > jgr->rid ) {
	++lastFirst;	
	continue;
      }
      else
	abort(); 
    }
    return nvisited;
    
    //this means end of file
    //bcf_destroy(v);
}

/**
 * Flush records.
 */
bcf1_t* JointGenotypeBlockReader::flush_variant(int32_t variantIndex, bcf_hdr_t* hdr, sex_ploidy_map& spmap) {
  // if it is the start of the block, read the entire files into PL/AD space
  int32_t nsamples = (int32_t) sample_names.size();
  //error("nsamples = %d", nsamples);
  if ( variantIndex == blockStarts[blockIDs[variantIndex]] ) {
    if ( variantIndex == 0 ) {
      blockPLs = (uint8_t*) malloc( sizeof(uint8_t) * unit * nsamples * 3 );
      blockADs = (uint8_t*) malloc( sizeof(uint8_t) * unit * nsamples * 3 );

      if ( ( blockPLs == NULL ) || ( blockADs == NULL ) )
	error("Cannot create a consecutive memory for storing %d x %d genotypes", unit, nsamples);
    }
    int32_t bid = blockIDs[variantIndex];
    int64_t bsize = blockEnds[bid] - blockStarts[bid] + 1;

    //htsFormat fmt;
    //fmt.category = unknown_category;
    //fmt.format = binary_format;
    //fmt.compression = no_compression;
    
    //htsFile* fp = hts_open_format(blockFNs[bid].c_str(), "r", &fmt);
    //htsFile* fp = hts_open(blockFNs[bid].c_str(), "rb");
    gzFile fp = gzopen(blockFNs[bid].c_str(), "rb");
    if ( fp == NULL )
      error("Cannot open the temporary file %s for reading", blockFNs[bid].c_str());

    int32_t sum_n_reads = 0;
    //notice("Flushing variant %d");
    int64_t nread = 0;
    for(int64_t s=0; s < nsamples; ++s) {
      for(int64_t i=0; i < bsize; ++i) {
	//notice("Flushing variant %d, s=%d, i=%d", variantIndex, s, i);
	//if ( bgzf_read(fp->fp.bgzf, blockPLs + (s + i * nsamples) * 3, sizeof(uint8_t)*3) < 0 )
	if ( (nread = gzread(fp, blockPLs + (s + i * nsamples) * 3, sizeof(uint8_t)*3)) != 3 )
	  error("Cannot read from temporary file %s, sample %d, variant %d/PL, returning %d", blockFNs[bid].c_str(), s, i, nread);
	sum_n_reads += nread;
	//if ( bgzf_read(fp->fp.bgzf, blockADs + (s + i * nsamples) * 3, sizeof(uint8_t)*3) < 0 )
	if ( (nread = gzread(fp, blockADs + (s + i * nsamples) * 3, sizeof(uint8_t)*3)) != 3 ) {
	  if ( ( gzeof(fp) == 0 ) || ( nread != 0 ) )
	    error("Cannot read from temporary file %s, sample %d, variant %d/AD, returning %d", blockFNs[bid].c_str(), s, i, nread);
	}
	sum_n_reads += nread;
      }
    }

    notice("Block %d - Size - %d, Reads %d", bid, bsize, sum_n_reads);

    char tmpchar;
    if ( gzread(fp, &tmpchar, sizeof(uint8_t)) != 0 )
      error("EOF has not been reached for temporary file %s. n_reads = %d", blockFNs[bid].c_str(), sum_n_reads);

    //if ( hts_close(fp) )
    if ( gzclose(fp) )
      error("Cannot close the temporary file %s after reading", blockFNs[bid].c_str());
  }
  //notice("Flushing variant %d", variantIndex);
  int64_t offset = (int64_t) (variantIndex - blockStarts[blockIDs[variantIndex]]) * (int64_t) nsamples * 3;
  if ( pFreqEst == NULL ) {
    // TODO: normalize the eigenvectors
    //pSVD = new Eigen::BDCSVD<Eigen::MatrixXd>(eV, Eigen::ComputeThinU | Eigen::ComputeThinV);
    pFreqEst = new frequency_estimator(&eV);
  }
  //frequency_estimator freqEst(pSVD, 1.0);
  return gRecords[variantIndex]->flush_variant(hdr, spmap, blockPLs + offset, blockADs + offset, pFreqEst);
}

void JointGenotypeBlockReader::write_header(BCFOrderedWriter* odw, bool printTmpInfo) {
  // contig and sample names must be added beforehand
  // write eigenvectors here
  bcf_hdr_append(odw->hdr, "##INFO=<ID=AVGDP,Number=1,Type=Float,Description=\"Average Depth per Sample\">\n");	
  bcf_hdr_append(odw->hdr, "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Alternate Allele Counts\">\n");
  bcf_hdr_append(odw->hdr, "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total Number Allele Counts\">\n");
  //bcf_hdr_append(odw->hdr, "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Reads\">\n");
  bcf_hdr_append(odw->hdr, "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Alternate Allele Frequency from Best-guess Genotypes\">\n");
  bcf_hdr_append(odw->hdr, "##INFO=<ID=GC,Number=G,Type=Integer,Description=\"Genotype Counts\">\n");
  bcf_hdr_append(odw->hdr, "##INFO=<ID=GN,Number=1,Type=Integer,Description=\"Total Number of Genotypes\">\n");
  bcf_hdr_append(odw->hdr, "##INFO=<ID=HWEAF_P,Number=A,Type=Float,Description=\"Genotype likelihood based pooled allele frequency assuming HWE\">\n");
  bcf_hdr_append(odw->hdr, "##INFO=<ID=FIBC_P,Number=1,Type=Float,Description=\"Pooled F-statistic (inbreeding coefficient) calculated from genotype likelihoods\">\n");	
  bcf_hdr_append(odw->hdr, "##INFO=<ID=HWE_SLP_P,Number=1,Type=Float,Description=\"Signed log p-values testing statistics based Hardy Weinberg Equilibrium Test ignoring population substructure\">\n");
  bcf_hdr_append(odw->hdr, "##INFO=<ID=FIBC_I,Number=1,Type=Float,Description=\"Popultion-structure-adjusted F-statistic (inbreeding coefficient) calculated from genotype likelihoods\">\n");	
  bcf_hdr_append(odw->hdr, "##INFO=<ID=HWE_SLP_I,Number=1,Type=Float,Description=\"Signed log p-values testing statistics based Hardy Weinberg Equilibrium Test adjusting for population substructure\">\n");
  bcf_hdr_append(odw->hdr, "##INFO=<ID=MAX_IF,Number=1,Type=Float,Description=\"Maximum individual-specific allele frequencies\">\n");
  bcf_hdr_append(odw->hdr, "##INFO=<ID=MIN_IF,Number=1,Type=Float,Description=\"Minimum individual-specific allele frequencies\">\n");
  bcf_hdr_append(odw->hdr, "##INFO=<ID=BETA_IF,Number=%d,Type=Float,Description=\"Coefficients for intercept and each eigenvector to obtain ISAF\">\n");
  bcf_hdr_append(odw->hdr, "##INFO=<ID=ABE,Number=1,Type=Float,Description=\"Expected allele Balance towards Reference Allele on Heterozygous Sites\">\n");
  bcf_hdr_append(odw->hdr, "##INFO=<ID=ABZ,Number=1,Type=Float,Description=\"Average Z-scores of Allele Balance towards Reference Allele on Heterozygous Sites\">\n");	
  bcf_hdr_append(odw->hdr, "##INFO=<ID=NS_NREF,Number=1,Type=Integer,Description=\"Number of samples with non-reference reads\">\n");
  bcf_hdr_append(odw->hdr, "##INFO=<ID=BQZ,Number=1,Type=Float,Description=\"Correlation between base quality and alleles\">\n");
  //bcf_hdr_append(odw->hdr, "##INFO=<ID=MQZ,Number=1,Type=Float,Description=\"Correlation between mapping quality and alleles\">\n");
  bcf_hdr_append(odw->hdr, "##INFO=<ID=CYZ,Number=1,Type=Float,Description=\"Correlation between cycle and alleles\">\n");
  bcf_hdr_append(odw->hdr, "##INFO=<ID=STZ,Number=1,Type=Float,Description=\"Correlation between strand and alleles\">\n");
  bcf_hdr_append(odw->hdr, "##INFO=<ID=NMZ,Number=1,Type=Float,Description=\"Correlation between mismatch counts per read and alleles\">\n");
  bcf_hdr_append(odw->hdr, "##INFO=<ID=IOR,Number=1,Type=Float,Description=\"Inflated rate of observing of other alleles in log10 scale\">\n");
  bcf_hdr_append(odw->hdr, "##INFO=<ID=NM0,Number=1,Type=Float,Description=\"Average number of mismatches in the reads with ref alleles\">\n");	
  bcf_hdr_append(odw->hdr, "##INFO=<ID=NM1,Number=1,Type=Float,Description=\"Average number of mismatches in the reads with non-ref alleles\">\n");
  if ( printTmpInfo )
    bcf_hdr_append(odw->hdr, "##INFO=<ID=FLT20,Number=20,Type=Float,Description=\"20 sufficient statistics for enabling hierarchical calls\">\n");  
  bcf_hdr_append(odw->hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
  bcf_hdr_append(odw->hdr, "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n");
  bcf_hdr_append(odw->hdr, "##FORMAT=<ID=AD,Number=A,Type=Integer,Description=\"Allele Depth\">\n");
  bcf_hdr_append(odw->hdr, "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total Allele Depth, including unidentifiable alleles\">\n");
  bcf_hdr_append(odw->hdr, "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Phred-scale Genotype Likelihoods\">\n");

  bcf_hdr_append(odw->hdr, "##FILTER=<ID=overlap_snp,Description=\"Overlaps with snp\">");
  bcf_hdr_append(odw->hdr, "##FILTER=<ID=overlap_indel,Description=\"Overlaps with indel\">");
  bcf_hdr_append(odw->hdr, "##FILTER=<ID=overlap_vntr,Description=\"Overlaps with VNTR\">");

  odw->write_hdr();
}
