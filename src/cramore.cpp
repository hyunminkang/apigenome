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
int32_t cmdCramSparseGenotype(int32_t argc, char** argv);
int32_t cmdCramCustomSort(int32_t argc, char** argv);
int32_t cmdCramFlagStat(int32_t argc, char** argv);
int32_t cmdCramCompareBQs(int32_t argc, char** argv);

int32_t cmdScMultinomEM(int32_t argc, char** argv);
int32_t cmdScMapSTAMPs(int32_t argc, char** argv);
int32_t cmdScKallistoCount(int32_t argc, char** argv);

int32_t cmdVcfMendelDupConc(int32_t argc, char** argv);
int32_t cmdVcfSampleSummary(int32_t argc, char** argv);
int32_t cmdVcfSqueeze(int32_t argc, char** argv);

int32_t cmdCramProcapDetect(int32_t argc, char** argv);


/*
typedef struct {
  int min_baseQ, tid, max_bases;
  uint16_t * bases;
  samFile* fp;
  bam_hdr_t* h;
  char* ref;
  int len;
  faidx_t *fai;
  errmod_t *em;
} ct_t;

static int read_aln(void* data, bam1_t* b) {
  int ret;
  ct_t *g = (ct_t*) data;
  while(1) {
    ret = sam_read(g->fp, g->h, b);
    if ( ret < 0 ) break;
    if ( b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP) ) continue;
    if ( g->fai && b->core.tid >= 0 ) {
      if (b->core.tid != g->tid) { // then load the sequence
	free(g->ref);
	g->ref = fai_fetch(g->fai, g->h->target_name[b->core.tid], &g->len);
	g->tid = b->core.tid;
      }
      bam_prob_realn_core(b, g->ref, g->len, 1<<1|1);
    }
    break;
  }
  return ret;
}

int32_t runGeneForest(int32_t argc, char** argv) {
  std::string inVcf;
  std::string mapf;
  std::string reg;
  std::string outf;

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input Sites", NULL)
    LONG_STRING_PARAM("in-vcf",&inVcf, "Input VCF/BCF file")
    LONG_STRING_PARAM("map",&mapf, "Map file containing population of each individual")    
    LONG_STRING_PARAM("region",&reg,"Genomic region to focus on")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &out, "Output VCF file name")
    LONG_INT_PARAM("verbose",&verbose,"Frequency of verbose output (1/n)")    
  END_LONG_PARAMS();  

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();


  // sanity check of input arguments
  if ( inVcf.empty() || mapf.empty() || out.empty() ) {
    error("--in-vcf, --map, --out are required parameters");
  }
  
  htsFile* fp = hts_open(mapf.c_str(),"r");
  if ( fp == NULL )
    error("Cannot open file %s for reading", mapf.c_str());

  kstring_t str = {0,0,0};
  int32_t lstr = 0;
  int32_t* flds = NULL;
  int32_t nflds = 0;

  std::map<std::string, std::string> id2pop;
  std::map<std::string, std::string> id2con;  

  notice("Started Reading sample map information");
  
  while( ( lstr = hts_getline(fp, KS_SEP_LINE, &str) ) >= 0 ) {
    flds = ksplit(&str, 0, &nflds); // expects ID, POPULATION, CONTINENT
    if ( nflds != 3 )
      error("Expected 3 columns from %s but observed %d", mapf.c_str(), nflds);

    id2pop[flds[0]] = flds[1];
    id2con[flds[0]] = flds[2];
  }

  notice("Started Reading site information from VCF file");

  std::vector<GenomeInterval> intervals;
  if ( !reg.empty() ) {
    parse_intervals(intervals, "", reg);
  }
  BCFOrderedReader odr(inVcf, intervals);
  bcf1_t* iv = bcf_init();
  
  uint32_t *fls = NULL, *frs = NULL, *rls = NULL, *rrs = NULL, *als = NULL;
  uint32_t nfl = 0, nfr = 0, nrl = 0, nrr = 0, nal = 0;

  int32_t nhaps = bcf_hdr_nsamples(odr.hdr);

  std::vector<int32_t> lois(nhaps,0);
  std::vector<int32_t> rois(nhaps,0);  
    
  while( odr.read(iv) ) {
    bcf_unpack(iv, BCF_UN_ALL);

    if ( bcf_get_format_int32(odr.hdr, iv, "RR", &rrs, &nrr) < 0 )
      error("Cannot parse RR field");

    if ( bcf_get_format_int32(odr.hdr, iv, "RL", &rrs, &nrr) < 0 )
      error("Cannot parse RL field");

    if ( bcf_get_format_int32(odr.hdr, iv, "FR", &rrs, &nrr) < 0 )
      error("Cannot parse FR field");

    if ( bcf_get_format_int32(odr.hdr, iv, "FL", &rrs, &nrr) < 0 )
      error("Cannot parse FL field");

    if ( bcf_get_format_int32(odr.hdr, iv, "AL", &rrs, &nrr) < 0 )
      error("Cannot parse AL field");

    // first, order the individuals based on ranks - should be linear time
    for(int32_t i=0; i < nhaps; ++i) {
      rois[rrs[i]] = i;
      lois[lrs[i]] = i;
    }

    // next, order by matching length - should be nlog(n)
  }
}
*/

int32_t main(int32_t argc, char** argv) {
  if ( argc < 2 ) {
    printf("cramore -- Fast analytic tools for analyzing and manipulating SAM/BAM/CRAM files\n");
    printf("Copyright (c) 2016 Hyun Min Kang and Adrian Tan\n");
    printf("Usage : %s [command] [options]\n",argv[0]);
    printf("\tType one of the following commands below to get detailed usage. Type %s [command] -help for more detailed descriptions\n",argv[0]);
    printf("\t%s cram-sparse-genotypes [options] : Sparse genotyping from BAM files\n",argv[0]);
    printf("\t%s kallisto-count [options] : Count digital expressions from kallisto-aligned dropseq data\n",argv[0]);    
    printf("\t%s sc-map-stamps [options] : Map STAMPs in DropSeq cell-UMI barcode FASTQ files\n",argv[0]);
    printf("\t%s sc-multinom-em [options] : Perform multinomial EM algorithm for digital expression data\n",argv[0]);        
    printf("\t%s vcf-mendel-dup-conc [options] : Mendelian/Duplicate genotype concordance\n",argv[0]);
    printf("\t%s vcf-sample-summary [options] : Sample level summary of VCF file\n",argv[0]);
    printf("\t%s vcf-squeeze [options] : Simplify VCF file by stripped out unessential fields\n",argv[0]);
    printf("\t%s verify-pair-id [options] : Verify identify of cell barcodes, including doublets from SAM/BAM/CRAM files\n",argv[0]);        
  }
  else {
    std::string cmd(argv[1]);
    if ( ( cmd == "kallisto-count" ) || ( cmd == "sc-kallisto-count" ) ) {
      return cmdScKallistoCount(argc-1,argv+1);
    }
    else if ( cmd == "sc-map-stamps" ) {
      return cmdScMapSTAMPs(argc-1,argv+1);
    }    
    else if ( ( cmd == "verify-pair-id" ) || ( cmd == "cram-verify-pair-id" ) ) {
      return cmdCramVerifyPairID(argc-1,argv+1);  
    }
    else if ( cmd == "cram-sparse-genotype" ) {
      return cmdCramSparseGenotype(argc-1,argv+1);
    }
    else if ( cmd == "vcf-mendel-dup-conc") {
      return cmdVcfMendelDupConc(argc-1,argv+1);      
    }
    else if ( cmd == "vcf-sample-summary") {
      return cmdVcfSampleSummary(argc-1,argv+1);      
    }
    else if ( cmd == "vcf-squeeze") {
      return cmdVcfSqueeze(argc-1,argv+1);      
    }
    else if ( ( cmd == "sc-multinom-em") || ( cmd == "multinom-em") ) {
      return cmdScMultinomEM(argc-1, argv+1);
    }
    else if ( cmd == "cram-simul-contam") {
      return cmdCramSimulContam(argc-1, argv+1);
    }
    else if ( cmd == "cram-flagstat-all") {
      return cmdCramFlagStat(argc-1, argv+1);
    }
    else if ( cmd == "cram-procap-detect") {
      return cmdCramProcapDetect(argc-1, argv+1);      
    }
    else if ( cmd == "bwa-pipe" ) {
      //return runBwaPipe(argc-1, argv+1);
    }
    else {
      error("Unrecognized command %s\n",argv[1]);
    }
  }
  return 0;
}
