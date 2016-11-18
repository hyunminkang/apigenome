#include "cramore.h"

int32_t cmdCramFlagStat(int32_t argc, char** argv) {
  std::string inSam;   // SAM, BAM, or CRAM
  std::string outFile; // Output file name 

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Options for input SAM/BAM/CRAM", NULL)
    LONG_STRING_PARAM("sam",&inSam, "Input SAM/BAM/CRAM file. Must be sorted by coordinates and indexed")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out",&outFile,"Output file name")
  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  // load VCF files. This VCF should only contain hard genotypes in GT field
  samFile* in = NULL;
  bam_hdr_t *header = NULL;

  if ( inSam.empty() )
    error("[E:%s:%d %s] --sam parameter is missing",__FILE__,__LINE__,__FUNCTION__);  

  if ( ( in = sam_open(inSam.c_str(), "r") ) == 0 ) {
    error("[E:%s:%d %s] Cannot open file %s\n",__FILE__,__LINE__,__FUNCTION__,inSam.c_str());    
  }

  if ( ( header = sam_hdr_read(in) ) == 0 ) {
    error("[E:%s:%d %s] Cannot open header from %s\n",__FILE__,__LINE__,__FUNCTION__,inSam.c_str());
  }

  if ( outFile.empty() )
    error("[E:%s:%d %s] --out parameter is missing",__FILE__,__LINE__,__FUNCTION__);

  bam1_t *b = bam_init1();
  //int32_t r;    

  //kstring_t readseq = {0,0,0};
  //kstring_t readqual = {0,0,0};

  notice("Reading SAM/BAM/CRAM records\n");

  int32_t nFlags = 12;  
  int32_t nFlagCodes = 65536;
  uint64_t* flagCounts = (uint64_t*)calloc(nFlagCodes, sizeof(uint64_t));
  int32_t ret;
  int32_t nReads = 0;
  while( ( ret = sam_read1(in, header, b) ) >= 0 ) {
    ++(flagCounts[b->core.flag]);
    ++nReads;
  }

  notice("Finished reading %d SAM/BAM/CRAM records. Now writing the full flagstat map", nReads);

  const char* flagDesc[] = {"PAIRED_OR_MULTI","PROP_PAIR","THIS_UNMAPPED","NEXT_UNMAPPED","THIS_REVCOMP","NEXT_REVCOMP","FIRST_SEGMENT","LAST_SEGMENT","SECONDARY","QC_FAILED","DUPLICATE","SUPPLEMENTARY"};
  int32_t i, j, k;

  std::vector< int32_t > flagStatus(nFlags, 0); // 01 : has zero, 10 : has one, 11 : has both
  for(i=0; i < nFlagCodes; ++i) {
    if ( flagCounts[i] > 0 ) {
      for(j=0; j < nFlags; ++j) {
	bool isFlagJOn = ( ( i & ( 1 << j ) ) ? true : false );
	if ( isFlagJOn ) {
	  flagStatus[j] |= 0x02;
	}
	else {
	  flagStatus[j] |= 0x01;
	}
      }
    }    
  }

  std::map< int32_t, std::map<int32_t, uint64_t> > maskMap;
  
  std::vector< int32_t > possibleMasks(1, 0);
  for(i=0; i < nFlags; ++i) {
    if ( flagStatus[i] == 0x03 ) {  // both values exist
      int size = (int)possibleMasks.size();
      for(j=0; j < size; ++j) {
	possibleMasks.push_back(possibleMasks[j] | (1 << i));
      }
    }
    else if ( ( flagStatus[i] <= 0 ) || ( flagStatus[i] > 3 ) ) {
      error("[E:%s:%d %s] Cannot recognize flagStatus %d",__FILE__,__LINE__,__FUNCTION__,flagStatus[i]);
    }
  }

  notice("# Possible masks = %lu", possibleMasks.size());

  for(i=0; i < nFlagCodes; ++i) {
    if ( flagCounts[i] > 0 ) {
      for(j=0; j < (int)possibleMasks.size(); ++j) {
	k = (possibleMasks[j] | i);
	maskMap[possibleMasks[j]][k] += flagCounts[i];
      }
    }
  }
    
  htsFile* wout = hts_open(outFile.c_str(), "w");

  // print header info
  for(j=0; j < nFlags; ++j) 
    hprintf(wout, "##%03x\t%s\n", 1 << j, flagDesc[j]);
  hprintf(wout, "#READ_COUNTS\tMASK\tMASKED\n");
  for(std::map< int32_t, std::map<int32_t, uint64_t> >::iterator it = maskMap.begin(); it != maskMap.end(); ++it) {
    for( std::map<int32_t, uint64_t>::iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
      hprintf(wout, "%15llu\t%03x\t%03x\n", it2->second, it->first, it2->first);
    }
    hprintf(wout,"\n");
  }

  notice("Finished writing output files");

  hts_close(wout);
  free(flagCounts);
  
  return 0;
}
