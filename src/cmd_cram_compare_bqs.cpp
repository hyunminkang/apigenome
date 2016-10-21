#include "cramore.h"

int32_t cmdCramCompareBQs(int32_t argc, char** argv) {
  std::string inSam1;    // SAM, BAM, or CRAM
  std::string inSam2;    // SAM, BAM, or CRAM  
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
    LONG_MULTI_STRING_PARAM("excl-vcf",&exclVcfs, "Input VCF files")

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
  bcf1_t* iv = bcf_init();

  samFile* in = NULL;
  bam_hdr_t *header = NULL;

  if ( ( in = sam_open(inSam.c_str(), "r") ) == 0 ) {
    error("Cannot open file %s\n",inSam.c_str());    
  }

  if ( ( header = sam_hdr_read(in) ) == 0 ) {
    error("Cannot open header from %s\n",inSam.c_str());
  }

  if ( outPrefix.empty() )
    error("--out parameter is missing");

  bam1_t *b = bam_init1();
  //int32_t r;    

  int32_t nsamples = bcf_hdr_nsamples(odr.hdr);
  std::vector< std::vector<uint8_t> > v_gts;
  std::vector<double> v_afs;
  std::vector<int32_t> v_rids;
  std::vector<int32_t> v_poss;
  std::vector<char> v_refs;
  std::vector<char> v_alts;

  
  
  samFile* in = NULL;
  bam_hdr_t *header = NULL;

  if ( inSam.empty() )
    error("--sam parameter is missing");  

  if ( ( in = sam_open(inSam.c_str(), "r") ) == 0 ) {
    error("Cannot open file %s\n",inSam.c_str());    
  }

  if ( ( header = sam_hdr_read(in) ) == 0 ) {
    error("Cannot open header from %s\n",inSam.c_str());
  }

  if ( outFile.empty() )
    error("--out parameter is missing");

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
      error("Cannot recognize flagStatus %d",flagStatus[i]);
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
