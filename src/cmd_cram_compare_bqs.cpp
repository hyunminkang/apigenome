#include "cramore.h"

int32_t cmdCramCompareBQs(int32_t argc, char** argv) {
  std::string inSam1;    // SAM, BAM, or CRAM
  std::string inSam2;    // SAM, BAM, or CRAM
  std::string refFasta;  // reference fasta file
  std::string outPrefix; // Output file name
  std::vector<std::string> exclVcfs; // positions to exclude as VCF files
  bool covCycle = false;
  bool covChr   = false;
  bool covDinuc = false;
  bool covRG = false;

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Options for input SAM/BAM/CRAM", NULL)
    LONG_STRING_PARAM("sam1",&inSam1, "Input SAM/BAM/CRAM file. Must be sorted by coordinates and indexed")
    LONG_STRING_PARAM("sam2",&inSam2, "Input SAM/BAM/CRAM file. Must be sorted by coordinates and indexed")
    LONG_STRING_PAKRAM("ref",&refFasta, "Reference FASTA file")
    LONG_MULTI_STRING_PARAM("excl-vcf",&exclVcfs, "Input VCF files to exclude for empirical quality calculation")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_PARAM("cov-cycle",&covCycle,"Use cycle as covariate")
    LONG_PARAM("cov-chr",&covChr,"Use chr as covariate")
    LONG_PARAM("cov-dinuc",&covDinuc,"Use dinucleotide as covariate")
    LONG_PARAM("cov-rg",&covCycle,"Use readgroup as covariate")    
    LONG_STRING_PARAM("out",&outPrefix,"Output file prefix")
  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  // load VCF files. This VCF should only contain hard genotypes in GT field
  std::vector<GenomeInterval> intervals;
  BCFOrderedReader odr(inVcf, intervals);

  samFile* in1 = NULL, in2 = NULL;
  bam_hdr_t *header1 = NULL, *header2 = NULL;

  if ( refFasta.empty() || outPrefix.empty() || inSam1.empty() || inSam2.empty() )
    error("One of the required parameters, --sam1, --sam2, --ref, or --out is missing\n");

  fastaMask fm(refFasta);
  for(int32_t i=0; i < exclVcfs.size(); ++i)
    fm.maskVcf(exclVcfs[i].c_str());

  if ( ( in1 = sam_open(inSam1.c_str(), "r") ) == 0 ) {
    error("Cannot open file %s\n",inSam1.c_str());    
  }

  if ( ( header1 = sam_hdr_read(in1) ) == 0 ) {
    error("Cannot open header from %s\n",inSam1.c_str());
  }

  if ( ( in2 = sam_open(inSam2.c_str(), "r") ) == 0 ) {
    error("Cannot open file %s\n",inSam2.c_str());    
  }

  if ( ( header2 = sam_hdr_read(in2) ) == 0 ) {
    error("Cannot open header from %s\n",inSam2.c_str());
  }

  bam1_t *b1 = bam_init1();
  bam1_t *b2 = bam_init1();  

  htsFile* wout = hts_open(outFile.c_str(), "w");

  int32_t ret1, ret2;
  kstring_t qual1 = {0,0,0};
  kstring_t qual2 = {0,0,0};
  kstring_t cstr = {0,0,0};  
  while( ( ret1 = sam_read1(in1, header1, b1) ) >= 0 ) {
    ret2 = sam_read2(in2, header2, b2);
    if ( ret1 != ret2 )
      error("[E:%s:%d %s] ret1 = %d != ret2 = %d",__FILE__,__LINE__,__FUNCTION__,ret1,ret2);

    bam1_core_t *c1 = b1->core;
    int32_t rlen = c1->l_qseq;
    uint32_t cpos = c1->pos;
    int32_t tid = c1->tid;

    bam1_core_t *c2 = b2->core;
    if ( ( rlen != c2->l_qseq ) && ( cpos != c2->pos ) )
      error("[E:%s:%d %s] SAM records are not identical between the files",__FILE__,__LINE__,__FUNCTION__);

    bam_get_qual_string(b1, &qual1);
    bam_get_qual_string(b2, &qual2);
    bam_get_seq_string(b1, &read1);
    bam_get_seq_string(b2, &read2);

    char rc;

    if ( c1->n_cigar ) {
      uint32_t* cigar = bam_get_cigar(b1);
      int32_t rpos = 0;
      int32_t m1, m2;
      for(uint32_t i=0; i < c1->n_cigar; ++i) {
	char op = bam_cigar_opchr(cigar[i]);
	cstr.l = 0;
	kputw(bam_cigar_oplen(cigar[i]), &cstr);
	char* stop;
	uint32_t len = strtol(str.s, &stop, 10);
	assert(stop);

	if ( op == 'M' ) { // cpos..cposL-1
	  for(uint32_t i=cpos; i < cpos+len; ++i) {
	    switch(rc = fm.getMask(tid, i)) {
	    case 'A': case 'C': case 'G': case 'T':
	      m1 = ( read1[rpos] == rc ? 0 : 1 );
	      m2 = ( read1[rpos] == rc ? 0 : 1 );	      
	      break;
	    case 'a': case 'c': case 'g': case 't':
	      m1 = ( read1[rpos] == toupper(rc) ? 0 : 1 );
	      m2 = ( read1[rpos] == toupper(rc) ? 0 : 1 );	      	      
	      break;
	    default:
	      m1 = m2 = 2;
	    }
	  }
	}
	else if ( ( op == 'D' ) || ( op == 'N' ) ) {
	}
	else if ( ( op == 'S' ) || ( op == 'I' ) ) {
	}
      }
    }
  }  

  hts_close(wout);
  
  
  return 0;
}
