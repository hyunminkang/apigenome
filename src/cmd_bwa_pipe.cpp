#include "cmd_bwa_pipe.h"
#include "Error.h"
#include "bam_ordered_reader.h"
#include "params.h"
#include "utils.h"

// extract
// <

int32_t cmdBwaPipe(int32_t argc, char** argv) {
  std::string in("-");
  std::string out("-");
  std::string ref;
  bool outSam = true;
  bool outBam = false;
  bool outCram = false;

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input/Output Files", NULL)
    LONG_STRING_PARAM("in", &in, "Input BAM/SAM/CRAM file (each pair of reads appears in consecutive lines)")
    LONG_STRING_PARAM("out", &out, "Output BAM/SAM/CRAM files (default in SAM format, use --out-bam or --out-cram to change)")
    LONG_STRING_PARAM("ref", &ref, "Reference sequence FASTA file")    

    LONG_PARAM_GROUP("Output file format", NULL)
    LONG_EXCLUSIVE_PARAM("out-sam", &outSam, "Output format is SAM")
    LONG_EXCLUSIVE_PARAM("out-bam", &outBam, "Output format is BAM")
    LONG_EXCLUSIVE_PARAM("out-cram", &outCram, "Output format is CRAM")
  END_LONG_PARAMS();  

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  samFile* fpIn = sam_open(in.c_str(),"r");
  if ( fpIn == NULL )
    error("Cannot open file %s for reading",in.c_str());

  bam_hdr_t* hdr = sam_hdr_read(fpIn);
  if ( hdr == NULL )
    error("Cannot open header from %s",in.c_str());

  std::vector<bam1_t*> readPools;
  reads.push_back(bam_init1());
  int32_t curIdx = 0;
  while( ( ret = sam_read1(fpIn, hdr, readPools[curIdx]) ) >= 0 ) { // read as long as the read names are the same
    if ( curIdx == 0 )
      continue;   // if there is only one read in the pool, read more
  }

  return 0;
  
}
