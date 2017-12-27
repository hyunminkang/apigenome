//#include <iostream>
//#include <string>
//#include <map>
//#include <vector>
//#include <ctime>
//#include <cstdlib>
//#ifdef _OPENMP
//  #include <omp.h>
//#else
//  #define omp_get_thread_num() 0
//#endif  

#include "cramore.h"
#include "commands.h"

//#include "htslib/sam.h"
//#include "htslib/faidx.h"
//#include "htslib/kseq.h"
//#include "htslib/khash.h"
//#include "estimator.h"
//#include "dropseq.h"
#include "utils.h"

KHASH_MAP_INIT_STR(vdict, bcf_idinfo_t)
typedef khash_t(vdict) vdict_t;

int32_t cmdCramSimulContam(int32_t argc, char** argv);
int32_t cmdCramVerifyPairID(int32_t argc, char** argv);
int32_t cmdCramDemuxlet(int32_t argc, char** argv);
int32_t cmdCramFreemux(int32_t argc, char** argv);
int32_t cmdCramMuxPileup(int32_t argc, char** argv);
int32_t cmdCramSimuxlet(int32_t argc, char** argv);
int32_t cmdCramSparseGenotype(int32_t argc, char** argv);
int32_t cmdCramDenseGenotype(int32_t argc, char** argv);
int32_t cmdCramCustomSort(int32_t argc, char** argv);
int32_t cmdCramFlagStat(int32_t argc, char** argv);
int32_t cmdCramCompareBQs(int32_t argc, char** argv);
int32_t cmdCramContextIndelAnalysis(int32_t argc, char** argv);

int32_t cmdScMultinomEM(int32_t argc, char** argv);
int32_t cmdScMultinomGibbs(int32_t argc, char** argv);
int32_t cmdScMapSTAMPs(int32_t argc, char** argv);
int32_t cmdScKallistoCount(int32_t argc, char** argv);

int32_t cmdVcfMendelDupConc(int32_t argc, char** argv);
int32_t cmdVcfSampleSummary(int32_t argc, char** argv);
int32_t cmdVcfSqueeze(int32_t argc, char** argv);
int32_t cmdVcfDeltaSVM(int32_t argc, char** argv);
int32_t cmdVcfExtract(int32_t argc, char** argv);
int32_t cmdVcfSVD(int32_t argc, char** argv);
int32_t cmdVcfInferAncestry(int32_t argc, char** argv);
int32_t cmdVcfInferISAF(int32_t argc, char** argv);
int32_t cmdVcfUpdateSites(int32_t argc, char** argv);
int32_t cmdVcfPasteCalls(int32_t argc, char** argv);
int32_t cmdVcfMergeCandidateVariants(int32_t argc, char** argv);

int32_t cmdBedDeltaSVMTrain(int32_t argc, char** argv);
int32_t cmdBedMatchedShuffle(int32_t argc, char** argv);
int32_t cmdBedShuffle(int32_t argc, char** argv);

int32_t cmdCramProcapDetect(int32_t argc, char** argv);
int32_t cmdFastaGCContent(int32_t argc, char** argv);
int32_t cmdVcfNormalizeDepth(int32_t argc, char** argv);

int32_t cmdBgenToVcf(int32_t argc, char** argv);
int32_t cmdCramVerifyBam(int32_t argc, char** argv);
int32_t cmdCramUpdateRG(int32_t argc, char** argv);


int32_t main(int32_t argc, char** argv) {
  commandList cl;

  BEGIN_LONG_COMMANDS(longCommandlines)
    LONG_COMMAND_GROUP("Processing BCF/VCF", NULL)
    LONG_COMMAND("mendel-dup-conc",&cmdVcfMendelDupConc, "Mendelian concordance analysis")
    LONG_COMMAND("vcf-sample-summary",&cmdVcfSampleSummary, "Sample-level summary from BCF/VCF")
    LONG_COMMAND("vcf-squeeze",&cmdVcfSqueeze, "Squeeze genotype fields from BCF/VCF")
    LONG_COMMAND("vcf-update-sites",&cmdVcfUpdateSites, "Update VCF site information")
    LONG_COMMAND("vcf-paste-calls",&cmdVcfPasteCalls, "Paste VCF calls produced by cram-dense-genotypes")
    LONG_COMMAND("vcf-merge-candidate-variants",&cmdVcfMergeCandidateVariants, "Merge candidate VCF sites produced by vt discover2")    

    LONG_COMMAND_GROUP("Processing SAM/BAM/CRAM", NULL)
    LONG_COMMAND("cram-simul-contam", &cmdCramSimulContam, "Simulate contaminated sequence reads")
    LONG_COMMAND("cram-context-indel-analysis", &cmdCramContextIndelAnalysis, "Context-specific analysis of indels")
    LONG_COMMAND("cram-flagstat-all", &cmdCramFlagStat, "Comprehensive flagstat analysis")
    LONG_COMMAND("cram-update-rg", &cmdCramUpdateRG, "Update readgroup information")

    LONG_COMMAND("cram-verify-bam", &cmdCramVerifyBam, "Ancestry-agnostic contamination estimation")
    
    LONG_COMMAND_GROUP("Population genetic analysis", NULL)
    LONG_COMMAND("vcf-svd",&cmdVcfSVD, "Perform SVD on BCF/VCF")
    LONG_COMMAND("vcf-infer-ancestry",&cmdVcfInferAncestry, "Infer genetic ancestry from VCF")
    LONG_COMMAND("vcf-infer-isaf",&cmdVcfInferISAF, "Infer individual-specific allele frequencies")
    LONG_COMMAND("vcf-normalize-depth", &cmdVcfNormalizeDepth, "Calculate GC-normalized depths for each chromosome")
    LONG_COMMAND("bgen2vcf",&cmdBgenToVcf, "Convert BGEN format to BCF/VCF")
    
    LONG_COMMAND_GROUP("Variant calling and genotyping", NULL)
    LONG_COMMAND("spare-genotype",&cmdCramSparseGenotype,"Sparse genotyping tool")
    LONG_COMMAND("dense-genotype",&cmdCramDenseGenotype, "Dense genotyping tool")
    

    LONG_COMMAND_GROUP("Functional analysis", NULL)    
    LONG_COMMAND("vcf-delta-svm",&cmdVcfDeltaSVM, "Delta SVM from BCF/VCF using lsgkm")
    LONG_COMMAND("bed-delta-svm-train",&cmdBedDeltaSVMTrain, "Train deltaSVM models")
    LONG_COMMAND("vcf-extract",&cmdVcfExtract, "Extract specific sites from BCF/VCF")
    LONG_COMMAND("cram-procap-detect", &cmdCramProcapDetect, "Detect TSS from PRO-cap data")
    
    LONG_COMMAND_GROUP("Single cell analysis", NULL)
    LONG_COMMAND("demuxlet", &cmdCramDemuxlet, "Deconvolute sample identify of droplet-based sc-RNAseq")
    //LONG_COMMAND("freemux", &cmdCramFreemux, "Genotype-free deconvolution of RNAseq")
    LONG_COMMAND("mux-pileup", &cmdCramMuxPileup, "Produce pileup of dsc-RNAseq")
    LONG_COMMAND("simuxlet",   &cmdCramSimuxlet,  "Simulate multiplexed dsc-RNAseq droplets")
    LONG_COMMAND("kallisto-count", &cmdScKallistoCount, "Produce digital expression matrix from kallisto-aligned sequence reads")
    LONG_COMMAND("sc-map-stamps", &cmdScMapSTAMPs, "Produce STAMP-map from DropSeq FASTQ files")
    LONG_COMMAND("sc-multinom-em", &cmdScMultinomEM, "Multinomial EM clustering of single cell types")
    LONG_COMMAND("sc-multinom-gibbs", &cmdScMultinomGibbs, "Multinomial Gibbs sampling of single cell types")    

    LONG_COMMAND_GROUP("Other tools", NULL)
    LONG_COMMAND("bed-matched-shuffle",&cmdBedMatchedShuffle, "Shuffle BED regions adjusting for GC contents and repeats")
    LONG_COMMAND("bed-shuffle",&cmdBedShuffle, "Shuffle BED regions randomly")        
    LONG_COMMAND("fasta-gc-content", &cmdFastaGCContent, "Create GC content profile")
  END_LONG_COMMANDS();

  cl.Add(new longCommands("Available Commands", longCommandlines));
  
  if ( argc < 2 ) {
    printf("[cramore] -- Fast analytic tools for analyzing high-throughput genomic data\n\n");
    fprintf(stderr, " Copyright (c) 2009-2017 by Hyun Min Kang and Adrian Tan\n");
    fprintf(stderr, " Licensed under the Apache License v2.0 http://www.apache.org/licenses/\n\n");    
    fprintf(stderr, "To run a specific command      : %s [command] [options]\n",argv[0]);
    fprintf(stderr, "For detailed instructions, run : %s --help\n",argv[0]);        
    cl.Status();
    return 1;
  }
  else {
    if ( strcmp(argv[1],"--help") == 0 ) {
      cl.HelpMessage();
    }
    else {
      return cl.Read(argc, argv);
    }
  }
  return 0;
}
