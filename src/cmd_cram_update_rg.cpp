#include "cramore.h"
#include "tsv_reader.h"

int32_t cmdCramUpdateRG(int32_t argc, char** argv) {
  std::string inSam;   // SAM, BAM, or CRAM
  std::string outSam;  // Output file name
  std::string rgMap;   // List of RGs to update
  std::string fasta;   // Name of reference FASTA file
  int32_t verbose = 10000000;
  
  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input Options", NULL)
    LONG_STRING_PARAM("sam",&inSam, "Input SAM/BAM/CRAM file. Must be sorted by coordinates and indexed")
    LONG_STRING_PARAM("rg-map",&rgMap, "Map of readgroups with [REPRESENTITAIVE_RG] [RG1] [RG2] ...")
    LONG_STRING_PARAM("ref", &fasta, "FASTA file (with .fai index) of the reference sequence")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out",&outSam,"Output file name")
    LONG_INT_PARAM("verbose",&verbose,"Verbose parameter")
  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  if ( inSam.empty() || outSam.empty() || rgMap.empty() || fasta.empty() )
    error("[E:%s:%d %s] required parameters --in, --out --rg-map, --ref is missing",__FILE__,__LINE__,__FUNCTION__);

  // load VCF files. This VCF should only contain hard genotypes in GT field
  samFile* in = NULL;
  samFile* out = NULL;
  bam_hdr_t *header = NULL;
  bam_hdr_t *outHeader = NULL;

  std::map<std::string,std::string> mRG;
  std::map<std::string,bool> rgs2remove;
  std::map<std::string,bool> rgs2keep;

  const char* fn_list = samfaipath(fasta.c_str());
  
  tsv_reader tsvf(rgMap.c_str());
  while( tsvf.read_line() > 0 ) {
    std::string keyRG(tsvf.str_field_at(0));
    rgs2keep.insert(std::pair<std::string,bool>(keyRG,false));
    mRG[keyRG] = keyRG;
    for(int32_t i=1; i < tsvf.nfields; ++i) {
      mRG[tsvf.str_field_at(i)] = keyRG;
      if ( keyRG.compare(tsvf.str_field_at(i)) != 0 )
	rgs2remove.insert(std::pair<std::string,bool>(tsvf.str_field_at(i),false));
    }
  }
    
  if ( inSam.empty() )
    error("[E:%s:%d %s] --sam parameter is missing",__FILE__,__LINE__,__FUNCTION__);  

  if ( ( in = sam_open(inSam.c_str(), "r") ) == 0 ) {
    error("[E:%s:%d %s] Cannot open file %s\n",__FILE__,__LINE__,__FUNCTION__,inSam.c_str());    
  }
  if ( hts_set_fai_filename(in, fn_list) != 0 ) 
    error("Failed to use reference %s", fn_list);  

  if ( ( header = sam_hdr_read(in) ) == 0 ) {
    error("[E:%s:%d %s] Cannot open header from %s\n",__FILE__,__LINE__,__FUNCTION__,inSam.c_str());
  }
  
  // make a new header
  const char* mode = ( outSam.substr( outSam.size() - 5 ).compare(".cram") == 0 ? "wc" : ( outSam.substr( outSam.size() - 5 ).compare(".bam") == 0 ? "wb" : "w" ) );
  if ( ( out = sam_open(outSam.c_str(), mode ) ) == 0 ) {
    error("[E:%s:%d %s] Cannot open file %s\n",__FILE__,__LINE__,__FUNCTION__,outSam.c_str());        
  }
  else {
    if ( hts_set_fai_filename(out, fn_list) != 0 ) 
      error("Failed to use reference %s", fn_list);
    
    kstring_t str = {0, 0, 0};
    char *cp, *line;
    int32_t j, l;

    cp = header->text;
    
    for(l = 0; l+3 < (int32_t)header->l_text; ++l) {
      line = &cp[l];
      if (!(cp[l] == '@' && cp[l+1] == 'R' && cp[l+2] == 'G')) { // if not starting with RG
	while (l < (int32_t)header->l_text && cp[l] != '\n') {
	  kputc(cp[l], &str);
	  l++;
	}
	if ( l < (int32_t)header->l_text ) kputc(cp[l], &str);  // copy each character one by one
      }
      else {  // if start with RG, find ID
	while (cp[l] != '\n') {
	  while (l < (int32_t)header->l_text && cp[l] != '\n' && cp[l] != '\t')
	    l++;
	  if (l+4 < (int32_t)header->l_text && cp[l+1] == 'I' && cp[l+2] == 'D')  // found ID:
	    break;
        }
        if (cp[l] == '\n')
	  error("Error in parsing header line %s", line);
		
        l = (j = l+4);
        while (l < (int32_t)header->l_text && cp[l] != '\n' && cp[l] != '\t')
	  l++;

        char *id = (char*)malloc(l-j+1);
        strncpy(id, &cp[j], l-j);
        id[l-j] = 0;
	
	// check whether it is a valid RG ID
	if ( rgs2remove.find(id) != rgs2remove.end() ) {
	  if ( rgs2remove[id] )
	    error("@RG header ID %s is seen multiple times",id); 
	  notice("Skipping to write @RG header ID %s", id);
	  rgs2remove[id] = true;

	  while (l < (int32_t)header->l_text && cp[l] != '\n') {
	    l++;
	  }
	}
	else if ( rgs2keep.find(id) != rgs2keep.end() ) {
	  if ( rgs2keep[id] )
	    error("@RG header ID %s is seen multiple times",id); 	    
	  notice("Writing @RG header ID %s", id);
	  rgs2keep[id] = true;
	  
	  // reset l again
	  l = line - cp;
	  while (l < (int32_t)header->l_text && cp[l] != '\n') {
	    kputc(cp[l], &str);
	    l++;
	  }
	  if ( l < (int32_t)header->l_text ) kputc(cp[l], &str);  // copy each character one by one
	}
	else {
	  error("@RG header ID %s is not listed in the file %s", id, rgMap.c_str());
	}

	free(id);
      }
    }

    for( std::map<std::string,bool>::iterator it = rgs2keep.begin(); it != rgs2keep.end(); ++it ) {
      if ( !it->second )
	error("Readgroup ID %s in %s was not observed in the header", it->first.c_str(), rgMap.c_str());
    }

    for( std::map<std::string,bool>::iterator it = rgs2remove.begin(); it != rgs2remove.end(); ++it ) {
      if ( !it->second )
	error("Readgroup ID %s in %s was not observed in the header", it->first.c_str(), rgMap.c_str());
    }

    //notice("Header length = %d", str.l);
    //notice("Header string = %s", str.s);    
    outHeader = sam_hdr_parse(str.l, str.s);
    outHeader->l_text = str.l; outHeader->text = str.s;
    outHeader = sam_hdr_sanitise(outHeader);
    
    if ( outHeader == NULL )
      error("Cannot parse header text %s",str.s);

    if ( sam_hdr_write(out, outHeader) < 0 ) {
      error("Cannot write header");      
    }
  }
  

  bam1_t *b = bam_init1();

  notice("Reading SAM/BAM/CRAM records");

  char rgtag[2] = {'R','G'};
  //int32_t ret;
  int32_t i;
  for( i=0; sam_read1(in, header, b) >= 0; ++i ) {
    if ( i % verbose == 0 )
      notice("Processing %d records at %s:%d", i, (b->core.flag & BAM_FUNMAP) ? "*" : bam_get_chrom(header,b), b->core.pos);
    // get RG ID
    uint8_t *uid = (uint8_t*) bam_aux_get(b, rgtag);
    const char* sid = ( ( uid != NULL ) && ( *uid == 'Z' ) ) ?  bam_aux2Z(uid) : NULL;
    if ( ( sid != NULL ) && ( rgs2remove.find(sid) != rgs2remove.end() ) ) { // need to update;
      bam_aux_update_str(b, rgtag, (int32_t)mRG[sid].size()+1, mRG[sid].c_str());
    }
    if ( sam_write1(out, outHeader, b) < 0 )
      error("Cannot write SAM record");
  }
  notice("Finished processing %d records");
  hts_close(out);
  hts_close(in);
  
  return 0;
}
