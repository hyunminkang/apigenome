#include "cramore.h"
#include "estimator.h"
#include "htslib/kseq.h"
#include "joint_genotype_block_reader.h"
#include "joint_genotype_block_record.h"
#include "bam_ordered_reader.h"

void bam_print_key_values(bam_hdr_t *h, bam1_t *s)
{
  const char* chrom = bam_get_chrom(h, s);
  uint32_t pos1 = bam_get_pos1(s);
  kstring_t seq = {0,0,0};
  bam_get_seq_string(s, &seq);
  uint32_t len = bam_get_l_qseq(s);
  kstring_t qual = {0,0,0};
  bam_get_qual_string(s, &qual);
  kstring_t cigar_string = {0,0,0};
  bam_get_cigar_string(s, &cigar_string);
  kstring_t cigar_expanded_string = {0,0,0};
  bam_get_cigar_expanded_string(s, &cigar_expanded_string);
  //uint16_t flag = bam_get_flag(s);
  uint32_t mapq = bam_get_mapq(s);
  
  uint8_t *aux;
  char* md = NULL;
  (aux=bam_aux_get(s, "MD")) &&  (md = bam_aux2Z(aux));
  
  std::cerr << "##################" << "\n";
  std::cerr << "chrom:pos: " << chrom << ":" << pos1 << "\n";
  std::cerr << "read     : " << seq.s << "\n";
  std::cerr << "qual     : " << qual.s << "\n";
  std::cerr << "cigar_str: " << cigar_string.s << "\n";
  std::cerr << "cigar    : " << cigar_expanded_string.s << "\n";
  std::cerr << "len      : " << len << "\n";
  std::cerr << "mapq     : " << mapq << "\n";
  std::cerr << "mpos1    : " << bam_get_mpos1(s) << "\n";
  std::cerr << "mtid     : " << bam_get_mtid(s) << "\n";
  std::cerr << "md       : " << (aux?md:"") << "\n";
  std::cerr << "##################" << "\n";
  
  if (seq.m) free(seq.s);
  if (qual.m) free(qual.s);
  if (cigar_string.m) free(cigar_string.s);
  if (cigar_expanded_string.m) free(cigar_expanded_string.s);
}

typedef struct {
  int32_t read_exclude_flag;
  int32_t read_mapq_cutoff;
  int32_t tid;
} filter_read_params_t;

static bool filter_read(bam_hdr_t* h, bam1_t *s, filter_read_params_t* param) {
  khiter_t k;
  int32_t ret;
  
  if(bam_get_flag(s) & param->read_exclude_flag) {
    //1. unmapped
    //2. secondary alignment
    //3. not passing QC
    //4. PCR or optical duplicate
    //++(param->no_exclude_flag_reads);
    return false;
  }
  
  if (bam_get_mapq(s) < param->read_mapq_cutoff) {
    //filter short aligments and those with too many indels (?)
    //++(param->no_low_mapq_reads);
    return false;
  }
  
  //*****************************************************************
  //should we have an assertion on the correctness of the bam record?
  //Is, Ds not sandwiched in M
  //leading and trailing Is - convert to S
  //no Ms!!!!!
  //*****************************************************************
  int32_t n_cigar_op = bam_get_n_cigar_op(s);
  if (n_cigar_op) {
    uint32_t *cigar = bam_get_cigar(s);
    bool seenM = false;
    int32_t last_opchr = '^';
    
    for (int32_t i = 0; i < n_cigar_op; ++i) {
      int32_t opchr = bam_cigar_opchr(cigar[i]);
      //int32_t oplen = bam_cigar_oplen(cigar[i]);
      if (opchr=='S') {
	if (i!=0 && i!=n_cigar_op-1) {
	  notice("S issue");
	  bam_print_key_values(h, s);
	  //++malformed_cigar;
	}
      }
      else if (opchr=='M') {
	seenM = true;
      }
      else if (opchr=='D') {
	if (last_opchr!='M' || (i<=n_cigar_op && bam_cigar_opchr(cigar[i+1])!='M')) {
	  notice("D issue");
	  //++no_malformed_del_cigars;
	  bam_print_key_values(h, s);
	}
      }
      else if (opchr=='I') {
	if (last_opchr!='M' || (i<n_cigar_op && bam_cigar_opchr(cigar[i+1])!='M')) {
	  if (last_opchr!='M') {
	    if (last_opchr!='^' && last_opchr!='S') {
	      notice("leading I issue\n");
	      bam_print_key_values(h, s);
	      //++no_malformed_ins_cigars;
	    }
	    else {
	      //++no_salvageable_ins_cigars;
	    }
	  }
	  else if (i==n_cigar_op-1) {
	    //++no_salvageable_ins_cigars;
	  }
	  else if (i==n_cigar_op-2 && (bam_cigar_opchr(cigar[i+1])=='S')) {
	    //++no_salvageable_ins_cigars;
	  }
	  else {
	    notice("trailing I issue");
	    bam_print_key_values(h, s);
	    //++no_malformed_ins_cigars;
	  }
	}
      }
      
      last_opchr = opchr;
    }
    
    if (!seenM) {
      notice("NO! M issue");
      bam_print_key_values(h, s);
      //++no_unaligned_cigars;
    }
  }
  
  return true;
}

int32_t cmdCramDenseGenotype(int32_t argc, char** argv) {
  std::string inVcf;
  std::vector<std::string> inCrams;
  std::string inCramList;
  std::string out;
  std::string output_tmp_prefix;
  std::string reg;
  //double gl_adj = 0.01;
  int32_t xStart = 2699520;
  int32_t xStop = 154931044;
  std::string smID;
  std::string xLabel("X");
  std::string yLabel("Y");
  std::string mtLabel("MT");
  std::string sexMap;
  int32_t unit = 1000000;
  int32_t capBQ = 40;
  //std::string refFasta;
  double minContam = 0.01;
  bool printTmpInfo = false;
  filter_read_params_t param;
  param.read_mapq_cutoff = 0;
  param.read_exclude_flag = 0x0704;
  param.tid = -1;

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input Sites", NULL)
    LONG_STRING_PARAM("in-vcf",&inVcf, "Input VCF/BCF file to genotype")
    LONG_STRING_PARAM("in-cram-list",&inCramList, "File containing input CRAM files in the order of [SM_ID] [CRAM_PATH] [CONTAM] [PC_1] [PC_2] ... [PC_N]")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &out, "Output VCF file name")
    LONG_STRING_PARAM("tmp-prefix", &output_tmp_prefix, "Prefix for temporary file (same to --out parameter by default)")    

    LONG_PARAM_GROUP("Key options", NULL)
    LONG_STRING_PARAM("region",&reg,"Genomic region to focus on")
    LONG_INT_PARAM("unit",&unit,"Maximum number of variants to stay in-memory")
    LONG_INT_PARAM("min-mq",&param.read_mapq_cutoff,"Minimum mapping quality threshold to consider")
    LONG_INT_PARAM("cap-bq",&capBQ,"Cap base quality over the threshold into a number")
    LONG_STRING_PARAM("sex-map",&sexMap, "Sex map file, containing ID and sex (1 for male and 2 for female) for each individual")
    LONG_DOUBLE_PARAM("min-contam",&minContam, "Minimum genotype likelihood adjustment factor at homozygous sites as Pr(1|0/0) or Pr(0|1/1)")    
    LONG_INT_PARAM("exclude-flag",&param.read_exclude_flag, "Flag to exclude reads")
    LONG_PARAM("print-tmp-info",&printTmpInfo,"Print temporary values INFO fields to allow merging")

    LONG_PARAM_GROUP("Sex Chromosomes",NULL)
    LONG_STRING_PARAM("xLabel", &xLabel, "Contig name for X chromosome")
    LONG_STRING_PARAM("yLabel", &yLabel, "Contig name for Y chromosome")
    LONG_STRING_PARAM("mtLabel", &mtLabel, "Contig name for MT chromosome")
    LONG_INT_PARAM("xStart", &xStart, "Start base position of non-PAR region in X chromosome")
    LONG_INT_PARAM("xStop",  &xStop,  "End base position of non-PAR region in X chromosome")    
  END_LONG_PARAMS();
  
  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();
  
  // sanity check of input arguments
  if ( inVcf.empty() || out.empty() || inCramList.empty() ) {
    error("[E:%s:%d %s] --in-vcf, --out, --in-cram are required parameters",__FILE__,__LINE__,__FUNCTION__);
  }

  if ( output_tmp_prefix.empty() ) output_tmp_prefix = out;

  int32_t nsamples = 0;
  int32_t ncols = 0;
  std::vector<std::string> sample_names;
  std::vector<std::string> cram_paths;
  std::vector<double> contams;
  std::vector< std::vector<double> > evecs;

  htsFile* fp = hts_open(inCramList.c_str(),"r");
  if ( fp == NULL )
    error("[E:%s:%d %s] Cannot open file %s for reading",__FILE__,__LINE__,__FUNCTION__, inCramList.c_str());
  
  int32_t lstr = 0;
  int32_t* fields = NULL;
  int32_t n = 0;
  kstring_t str = {0,0,0};      
  while( ( lstr = hts_getline(fp, KS_SEP_LINE, &str) ) >= 0 ) {
    fields = ksplit(&str, 0, &n);
    
    if ( ncols == 0 ) ncols = n;
    else if ( ncols != n ) error("The number of lines are inconsistent at line %d - %d before vs %d in this line", nsamples+1, ncols, n);
    
	
    if ( n < 2 )
      error("[E:%s:%d %s] in-cram-list file %s contains whitespace - # fields = %d, (%s, %s)",__FILE__,__LINE__,__FUNCTION__, inCramList.c_str(), n, str.s + fields[0], str.s + fields[1]);

    sample_names.push_back(std::string(str.s + fields[0]));
    cram_paths.push_back(std::string(str.s + fields[1]));
    if ( n >= 3 ) {
      contams.push_back(atof(str.s + fields[2]));
      if ( contams.back() < minContam ) contams.back() = minContam;
    }
    else {
      contams.push_back(minContam);
    }
    
    evecs.resize(nsamples+1);
    if ( n > 3 ) {
      evecs.back().resize(n-3);
      for(int32_t i=3; i < n; ++i) {
	evecs.back().at(i-3) = atof(str.s + fields[i]);
      }
    }
    ++nsamples;
  }
  hts_close(fp);

  std::vector<GenomeInterval> intervals;
  if ( !reg.empty() ) {
    parse_intervals(intervals, "", reg);
  }

  //std::vector<int32_t> vSex;
  //std::map<std::string,int> mSex;

  /*
  if ( !sexMap.empty() ) {
    htsFile *file = hts_open(sexMap.c_str(),"r");
    if ( file == NULL ) {
      fprintf(stderr,"ERROR: Cannot open %s\n",sexMap.c_str());
      exit(1);
    }
    kstring_t *s = &file->line;
    while( hts_getline(file,'\n',s) >= 0 ) {
      std::string ss = std::string(s->s);
      size_t idx = ss.find_first_of("\t ");
      if ( idx == std::string::npos ) {
	fprintf(stderr,"ERROR: Cannot parse line %s in %s\n",ss.c_str(), sexMap.c_str());
	exit(1);
      }
      std::string id = ss.substr(0, idx);
      int32_t sex = atoi(ss.substr(idx+1).c_str());
      
      if ( mSex.find(id) != mSex.end() ) {
	fprintf(stderr,"ERROR: Duplicate ID %s in %s\n",id.c_str(), sexMap.c_str());
	exit(1);	      
      }
      
      if ( sex == 0 ) {
	fprintf(stderr,"WARNING: Unknown sex for individual %s, assuming females\n",id.c_str());
	sex = 2;
      }
      else if ( sex > 2 ) {
	fprintf(stderr,"ERROR: Invalid sex %d for individual %s\n",sex,id.c_str());
	exit(1);
      }
      mSex[id] = sex;
    }
  } 
  */

  // load the VCF file first
  notice("Loading input VCF file %s in region %s", inVcf.c_str(), intervals[0].to_string().c_str());
  JointGenotypeBlockReader jgbr(inVcf, intervals, output_tmp_prefix, nsamples, unit, printTmpInfo);

  BCFOrderedWriter* odw = new BCFOrderedWriter(out);
  bcf_hdr_transfer_contigs(jgbr.odr->hdr, odw->hdr);

  //int32_t x_rid = bcf_hdr_name2id(jgbr.odr->hdr, xLabel.c_str());
  bam1_t* s = bam_init1();

  for(int32_t i=0; i < nsamples; ++i) {
    BAMOrderedReader odr(cram_paths[i], intervals);
    bam_hdr_t* h = odr.hdr;
    int64_t no_reads = 0;
    int64_t no_filt_reads = 0;

    if ( sample_names[i].compare(bam_hdr_get_sample_name(h)) != 0 )
      warning("The same name %s is different from the same name %s from %s. Continuing with the former one", sample_names[i].c_str(), bam_hdr_get_sample_name(h).c_str(), cram_paths[i].c_str());

    jgbr.set_sample(i, sample_names[i].c_str(), contams[i], evecs[i]);

    if ( jgbr.numVariants() > 0 ) {
      while( odr.read(s) ) {
	++no_reads;
	if ( !filter_read(h, s, &param) ) continue;

	++no_filt_reads;
	jgbr.process_read(h, s, i);
      }
    }

    if ( no_reads == 0 )
      warning("No read found in %d-th sample %s", i+1, sample_names[i].c_str());

    notice("Processed %d-th sample %s across %lld reads (%.1lf%% passed filter)", i+1, sample_names[i].c_str(), no_reads, no_filt_reads/(double)no_reads*100);
  
    jgbr.flush_sample(i);
    odr.close();
    //delete odr;
  }
  jgbr.close_blocks();

  bam_destroy1(s);

  for(int32_t i=0; i < nsamples; ++i) {
    bcf_hdr_add_sample(odw->hdr, sample_names[i].c_str());
  }
  bcf_hdr_add_sample(odw->hdr, NULL);

  int32_t nvariants = jgbr.numVariants();
  jgbr.write_header(odw, printTmpInfo);

  sex_ploidy_map spmap(xLabel, yLabel, mtLabel, xStart, xStop);
  spmap.load_sex_map_file(sexMap.empty() ? NULL : sexMap.c_str(), odw->hdr);

  for(int32_t i=0; i < nvariants; ++i) {
    if ( i % 1000 == 0 )
      notice("Writing %d variants to BCF/VCF file %s", i, out.c_str());
    bcf1_t* nv = jgbr.flush_variant(i, odw->hdr, spmap);
    odw->write(nv);
    bcf_destroy(nv);
  }
  odw->close();
  delete odw;

  notice("Finished writing %d variants to BCF/VCF file %s", nvariants, out.c_str());
  return 0;
}
