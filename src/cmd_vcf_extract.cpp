#include "cramore.h"
#include "bcf_filter_arg.h"
#include "bcf_variant_key.h"
#include "bcf_chunked_reader.h"
#include "bcf_ordered_writer.h"

int32_t cmdVcfExtract(int32_t argc, char** argv) {
  std::string siteVcf;
  std::string inVcf;
  std::string out;
  std::string refFasta;
  int32_t unit = INT_MAX;
  int32_t verbose = 10000;
  bool jumpFlag = false;

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input VCFs", NULL)
    LONG_STRING_PARAM("in",&inVcf, "Input VCF/BCF files")
    LONG_STRING_PARAM("site",&siteVcf, "VCF file containing the site to subset")
    LONG_STRING_PARAM("ref",&refFasta, "FASTA file containing the reference sequence")
    LONG_INT_PARAM("unit",&unit, "Unit for chunking the genome")
    LONG_PARAM("jump",&jumpFlag,"Jump using index")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &out, "Output VCF file name")
    LONG_INT_PARAM("verbose",&verbose,"Frequency of verbose output (1/n)")    
  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();
  
  // sanity check of input arguments
  if ( inVcf.empty() || out.empty() || siteVcf.empty() ) {
    error("[E:%s:%d %s] --in, --out are required parameters",__FILE__,__LINE__,__FUNCTION__);
  }

  //notice("foo");
  BCFChunkedReader cdrS(siteVcf.c_str());
  //notice("bar");  
  BCFChunkedReader cdrI(inVcf.c_str(), refFasta.c_str(), unit);
  //notice("goo");  

  std::set<variantKeyS> variants;
  std::vector<variantKeyS> variantList;
  bcf1_t* iv = bcf_init();
  for(int32_t i=0;cdrS.read(iv); ++i ) {
    if ( i % verbose == 0 )
      notice("Reading %d variants..", i+1);
    variantList.push_back(variantKeyS(cdrS.hdr, iv) );
    variants.insert(variantList.back());    
  }
  cdrS.close();

  notice("Finished loading %u variants to extract", variants.size());

  int32_t nVariants = (int32_t)variants.size();

  BCFOrderedWriter odw(out.c_str());
  odw.set_hdr(cdrI.hdr);
  odw.write_hdr();

  int32_t nRead = 0, nWritten = 0;  
  if ( jumpFlag ) {
    int32_t nvar = (int32_t)variantList.size();
    for(int32_t i=0; i < nvar; ++i) {
      if ( nRead % verbose == 0 )
	notice("Reading %d variants (%d intended) and writing %d variants at %s:%d", nRead, i, nWritten, variantList[i].chrom.c_str(), variantList[i].pos);
      
      if ( cdrI.jump_to(variantList[i].chrom.c_str(), variantList[i].pos) ) {
	//notice("foo");
	for(int32_t j=0; cdrI.read(iv); ++j) {
	  //notice("bar %d %d",iv->pos, variantList[i].pos);	  
	  ++nRead;
	  variantKeyS key(cdrI.hdr, iv);
	  if ( ( variantList[i].chrom.compare(bcf_get_chrom(cdrI.hdr,iv)) != 0 ) || ( iv->pos > variantList[i].pos ) )
	    break;

	  if ( variants.find(key) != variants.end() ) {
	    odw.write(iv);
	    ++nWritten;
	    variants.erase(key);	    
	  }
	}
      }
      else {
	warning("Cannot jump to %s:%d",variantList[i].chrom.c_str(), variantList[i].pos);
      }
    }
  }
  else {
    for( int32_t i=0; cdrI.read(iv); ++i ) {
      variantKeyS key(cdrI.hdr, iv);
      
      if ( i % verbose == 0 )
	notice("Reading %d variants and writing %d variants at %s:%d", nRead, nWritten, key.chrom.c_str(), key.pos);    
      
      if ( variants.find(key) != variants.end() ) {
	odw.write(iv);
	variants.erase(key);	
	++nWritten;
      }
      ++nRead;
    }
  }
  odw.close();
  cdrI.close();

  notice("Finishing writing %d variants, missing %d", nVariants, nVariants-nWritten);

  return 0;
}

