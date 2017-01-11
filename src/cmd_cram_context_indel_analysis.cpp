#include "cramore.h"

/*
class SamCigar {
public:
  std::vector<char> opchrs;
  std::vector<int32_t> oplens;

  samCigar(bam1_t* b) {
    int32_t n_cigar_op = bam_get_n_cigar_op(b);
    uint32_t* cigar = bam_get_cigar(b);
    opchrs.resize(n_cigar_op,'\0');
    oplens.resize(n_cigar_op,0);
    for(int32_t i=0; i < n_cigar_op; ++i) {
      opchrs[i] = bam_cigar_opchr(cigar[i]);
      oplens[i] = bam_cigar_oplen(cigar[i]);
    }
  }

  bool hasIndels() {
    for(int32_t i=0; i < opchrs.size(); ++i) {
      if ( ( opchrs[i] == 'I' ) || ( opchrs[i] == 'D' ) )
	return true;
    }
    return false;
  }
};

void displayAlignment(SamCigar& cigar, char* ref, char* seq, int32_t buffer, FILE* fp = stdout) {
  // display the reference genome sequence first
  int32_t lref = strlen(ref);
  int32_t lseq = strlen(seq);

  
  while( *ref ) 
    fputc(*(ref++), fp);
  fputc('\n');
  
}

*/

int32_t cmdCramContextIndelAnalysis(int32_t argc, char** argv) {
  std::string inSam; // SAM, BAM, or CRAM
  std::string region; // region to focus

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Options for input SAM/BAM/CRAM", NULL)
    LONG_STRING_PARAM("sam",&inSam, "Input SAM/BAM/CRAM file. Must be sorted by coordinates and indexed")
    LONG_STRING_PARAM("region",&region, "Region to focus on")    
  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  /*
  std::vector<GenomeInterval> intervals;    
  hts_idx_t *idx = sam_index_load(in, inSam.c_str());
  hts_itr_t *itr = bam_itr_querys(idx, header, reg.c_str());

  int32_t nReadsAll = 0, nReadsPass = 0;
  
  while( sam_itr_next(in, itr, b) >= 0 ) {
      ++nReadsAll;

      if ( b->core.flag & qcExclFlag ) continue;

      SamCigar samCigar(b);
      if ( samCigar.hasIndels ) {
	++nReadsPass;
      }
      else {
      }
  }
  sam_itr_destroy(itr);
  sam_idx_destroy(idx);
  */
  
  return 0;
}
