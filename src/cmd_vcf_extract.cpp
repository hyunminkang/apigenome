#include "cramore.h"
#include "bcf_filter_arg.h"
#include "bcf_variant_key.h"
#include "bcf_filtered_reader.h"
#include "bcf_ordered_reader.h"
#include "bcf_ordered_writer.h"
#include "tsv_reader.h"

int32_t cmdVcfExtract(int32_t argc, char** argv) {
  BCFFilteredReader bfr;
  std::string siteVcf;
  std::string out;
  std::string refFasta;
  std::string rsList;
  std::string varList;
  int32_t unit = INT_MAX;
  int32_t verbose = 10000;
  bool jumpFlag = false;
  bool rsSorted = false;

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input VCFs", NULL)
    LONG_STRING_PARAM("vcf",&bfr.bcf_file_name, "Input VCF/BCF files")
    LONG_STRING_PARAM("site",&siteVcf, "VCF file containing the site to subset")
    LONG_STRING_PARAM("rs-list",&rsList, "Text file containing the list of rsIDs to extract (must have --site as a dbSNP file)")

    LONG_PARAM_GROUP("Options to specify when chunking is used for input BCF/VCF", NULL)    
    LONG_STRING_PARAM("ref",&bfr.ref_file_name, "Reference FASTA file name (specify only when chunking is used)")
    LONG_INT_PARAM("unit",&bfr.unit, "Chunking unit in bp (specify only with --ref together")
    LONG_STRING_PARAM("interval",&bfr.interval_file_name, "Interval file name used for chunking (specify only when chunking is used without --ref")
    
    LONG_STRING_PARAM("region",&bfr.target_region, "Target region to focus on")      
    LONG_STRING_PARAM("var-list",&varList, "Text file containing the list of variant keys in [CHR]:[POS]:[REF]:[ALTs] format. --site input is not necessary")    
    LONG_STRING_PARAM("ref",&refFasta, "FASTA file containing the reference sequence")
    LONG_PARAM("rs-sorted",&rsSorted, "(--rs-list only) The site file is specialized (rsID-sorted) version, which does not follow a standard VCF format")    
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
  if ( out.empty() ) {
    error("[E:%s:%d %s] --out is a require parameters",__FILE__,__LINE__,__FUNCTION__);
  }

  if ( siteVcf.empty() && varList.empty() ) {
    error("[E:%s:%d %s] --site & --rs-list, or --var-list is a required parameter",__FILE__,__LINE__,__FUNCTION__);    
  }

  std::set<variantKeyS> variants;
  std::vector<variantKeyS> variantList;

  if ( !varList.empty() ) { // parse varList
    tsv_reader tsv_var_list(varList.c_str());
    while( tsv_var_list.read_line() > 0 ) {
      variantList.push_back(variantKeyS(tsv_var_list.str_field_at(0)) );
      variants.insert(variantList.back());
    }
    notice("Finished loading %u variants from %s", variants.size(), varList.c_str());        
  }

  std::set<int32_t> rsIDs; // list of rsIDs
  
  if ( !rsList.empty() ) { // parse rsIDs;
    tsv_reader tsv(rsList.c_str());
    while( tsv.read_line() > 0 ) {
      const char* s = tsv.str_field_at(0);
      if ( ( s[0] == 'r' ) && ( s[1] == 's' ) ) {
	rsIDs.insert(atoi(s+2));
      }
      else {
	rsIDs.insert(atoi(s));
      }
    }
    notice("Finished reading %u rsIDs", rsIDs.size());
  }


  if ( rsSorted ) { // search for sorted tabixed text file
    // Read the first record to find out the chromosome name

    notice("Reading rs-sorted dbSNP file %s", siteVcf.c_str());
    
    tsv_reader tsvS(siteVcf.c_str());
    int ncol = tsvS.read_line();
    if ( ncol == 0 )
      error("[E:%s] Cannot read rs-sorted site VCF %s", __PRETTY_FUNCTION__, siteVcf.c_str());
    std::string chr = tsvS.str_field_at(0);
    std::set<int32_t>::iterator it;
    for(it = rsIDs.begin(); it != rsIDs.end(); ++it) {
      //notice("Reading %s:%d", chr.c_str(), *it);
      
      if ( tsvS.jump_to(chr.c_str(), *it, *it) ) {
	//notice("Jumped to  %s:%d", chr.c_str(), *it);	
	if ( tsvS.read_line() > 6 ) {
	  // extract 0+2, 1+2, 3+2, 4+2 columns
	  if ( tsvS.int_field_at(1) != *it )
	    error("[E:%s] Queried rsID is %d but found %d", *it, tsvS.int_field_at(1));
	  //notice("Read line %d", tsvS.int_field_at(1));
	  
	  std::string key(tsvS.str_field_at(2));
	  key += ":";
	  key += tsvS.str_field_at(3);
	  key += ":";	  
	  key += tsvS.str_field_at(5);
	  key += ":";	  
	  key += tsvS.str_field_at(6);

	  //notice("Inserting %s", key.c_str());
	  
	  variantList.push_back(variantKeyS(key.c_str()));
	  variants.insert(variantList.back());

	  //notice("Finished inserting %s", key.c_str());
	}
	else
	  notice("rsID rs%d cannot be read from the rs-sorted dbSNP file", *it);
      }
      else
	notice("rsID rs%d cannot be found in the rs-sorted dbSNP file", *it);

      //notice("Inserted %d", tsvS.int_field_at(1));      
    }
    
    notice("Finished reading rs-sorted dbSNP file %s and found %u of %u", siteVcf.c_str(), variantList.size(), rsIDs.size());
  }
  else {
    std::vector<GenomeInterval> intervals;
    if ( !bfr.target_region.empty() ) {
      intervals.push_back( GenomeInterval(bfr.target_region) );
    }
    BCFOrderedReader odrS(siteVcf, intervals);
    bcf1_t* iv = bcf_init();

    for(int32_t i=0; odrS.read(iv) ; ++i ) {
      if ( i % verbose == 0 )
	notice("Reading %d variants..", i+1);

      //if ( bfrS.passed_vfilter() ) {
      //bcf1_t* iv = bfrS.cursor();
      if ( rsIDs.empty() || ( rsIDs.find(atoi(iv->d.id)) != rsIDs.end() ) ) {
	variantList.push_back(variantKeyS(odrS.hdr, iv) );
	variants.insert(variantList.back());
      }
	//}
    }
    bcf_destroy(iv);
    //bfrS.close();
    notice("Finished loading %u variants to extract", variants.size());    
  }

  //notice("bar");
  //BCFChunkedReader cdrI(inVcf.c_str(), refFasta.empty() ? NULL : refFasta.c_str(), unit);
  //notice("goo");
  bfr.init_params();

  int32_t nVariants = (int32_t)variants.size();

  BCFOrderedWriter odw(out.c_str());
  odw.set_hdr(bfr.cdr.hdr);
  odw.write_hdr();

  int32_t nRead = 0, nWritten = 0;  
  if ( jumpFlag ) {
    //int32_t nvar = variants.size(); //(int32_t)variantList.size();
    std::set<variantKeyS>::iterator it;
    //int32_t i=0;
    for( it = variants.begin(); it != variants.end(); ++it) {
      if ( nRead % verbose == 0 )
	notice("Reading %d/%u variants and writing %d variants at %s:%d", nRead, variants.size(), nWritten, it->chrom.c_str(), it->pos);
      
      if ( bfr.jump_to(it->chrom.c_str(), it->pos) ) {
	//notice("foo");
	for(int32_t j=0; bfr.read(); ++j) {
	  //bcf1_t* iv = bfr.cursor();	  
	  //notice("bar i=%d/%d, j=%d, nRead=%d, %d %d", i, nvar, j, nRead, iv->pos, variantList[i].pos);	  
	  ++nRead;
	  variantKeyS key(bfr.cdr.hdr, bfr.cursor());

	  while ( ( bfr.cursor()->pos > it->pos ) && ( it->chrom.compare(bcf_get_chrom(bfr.cdr.hdr, bfr.cursor())) == 0 ) ) {
	    ++it;
	  }

	  if ( it->chrom.compare(bcf_get_chrom(bfr.cdr.hdr, bfr.cursor())) != 0 )
	    break;
	  
	  
	  //if ( ( it->chrom.compare(bcf_get_chrom(bfr.cdr.hdr, bfr.cursor())) != 0 ) || ( bfr.cursor()->pos > it->pos ) )
	  //break;

	  //if ( variants.find(key) != variants.end() ) {
	  if ( key == *it ) {
	    odw.write(bfr.cursor());
	    ++nWritten;
	    //variants.erase(key);	    
	  }
	}
	//notice("goo");	
      }
      else {
	warning("Cannot jump to %s:%d",it->chrom.c_str(), it->pos);
      }
    }
    variants.clear();
  }
  else {
    for( int32_t i=0; bfr.read(); ++i ) {
      bcf1_t* iv = bfr.cursor();
      variantKeyS key(bfr.cdr.hdr, iv);
      
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
  bfr.cdr.close();

  notice("Finishing writing %d variants, missing %d", nWritten, nVariants-nWritten);

  return 0;
}

