#include "cramore.h"
#include "bcf_filtered_reader.h"
#include "bcf_ordered_writer.h"
//#include "ancestry_estimator.h"
#include <map>
#include <string>
#include <ctime>

int32_t cmdVcfUpdateSites(int32_t argc, char** argv) {
  BCFFilteredReader genoVcf;
  BCFFilteredReader siteVcf;
  BCFFilteredReader dbsnpVcf;
  std::vector<std::string> info2update;
  std::vector<std::string> info2remove;
  bool replaceFilter = false;
  bool replaceQual = false;
  bool replaceID = false;
  std::string outVcf;

  genoVcf.verbose = 10000;
  siteVcf.verbose = 100000;
  dbsnpVcf.verbose = INT_MAX;

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input Files", NULL)
    LONG_STRING_PARAM("vcf",&genoVcf.bcf_file_name, "Original VCF files to be updated")
    LONG_STRING_PARAM("site",&siteVcf.bcf_file_name, "Site VCF file containing desired information. The VCF files must have identical set of variants with file specified with --vcf option")
    LONG_STRING_PARAM("dbsnp",&dbsnpVcf.bcf_file_name, "dbSNP file containing the ID information")
    LONG_STRING_PARAM("region",&dbsnpVcf.target_region, "Target regions to focus on")

    LONG_PARAM_GROUP("Restrictions on Site VCFs", NULL)    
    LONG_STRING_PARAM("include-expr",&siteVcf.vfilt.include_expr, "Include sites for which expression is true")
    LONG_STRING_PARAM("exclude-expr",&siteVcf.vfilt.exclude_expr, "Exclude sites for which expression is true")    
    
    
    LONG_PARAM_GROUP("Options for modification", NULL)    
    LONG_MULTI_STRING_PARAM("md-info", &info2update, "Name of INFO field to add/modify")
    LONG_MULTI_STRING_PARAM("rm-info", &info2remove, "INFO field to add/update")    
    LONG_PARAM("replace-filter", &replaceFilter, "Replace the FILTER column with the new one")
    LONG_PARAM("replace-qual", &replaceQual, "Replace the QUAL column with the new one")
    LONG_PARAM("replace-id", &replaceQual, "Replace the ID column with the new one")     
    LONG_PARAM_GROUP("Output Files", NULL)
    LONG_STRING_PARAM("out",&outVcf, "Output VCF file")
    END_LONG_PARAMS();
  
  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();
  
  notice("Analysis Started");    
  
  // sanity check of input arguments
  if ( outVcf.empty() || genoVcf.bcf_file_name.empty() || siteVcf.bcf_file_name.empty() ) {
    error("[E:%s:%d %s] --geno, --site, --out are required parameters",__FILE__,__LINE__,__PRETTY_FUNCTION__);
  }

  bool hasdb = dbsnpVcf.bcf_file_name.empty() ? false : true;
  if ( hasdb ) {
    notice("Opening dbSNP file %s...", dbsnpVcf.bcf_file_name.c_str());
    //dbsnpVcf.target_region = siteVcf.target_region;
    dbsnpVcf.init_params();
  }

  siteVcf.init_params();
  genoVcf.init_params();

  BCFOrderedWriter odw(outVcf.c_str(), 0);
  odw.set_hdr(genoVcf.cdr.hdr);
  
  if ( replaceFilter ) {
    // Remove all FILTER fields from the original VCF
    //for(int32_t i=0; i < odw.hdr->nhrec; ++i) {
    for(int32_t i=odw.hdr->nhrec-1; i >=0; --i) {      
      bcf_hrec_t *hrec = odw.hdr->hrec[i];
      //notice("foo i=%d/%d, hrec=%x", i, odw.hdr->nhrec, hrec);
      if ( hrec->type == BCF_HL_FLT ) {
	int k = bcf_hrec_find_key(hrec,"ID");
	assert( k>=0 ); // this should always be true for valid VCFs
	//int hdr_id = bcf_hdr_id2int(odw.hdr, BCF_DT_ID, hrec->vals[k]);
	//notice("Removing %s", hrec->vals[k]);
	bcf_hdr_remove(odw.hdr, BCF_HL_FLT, hrec->vals[k]);
	//bcf_remove_filter(odw.hdr, gv, hdr_id, 0);
      }
    }

    // Add all filter columns from the new VCF
    for(int32_t i=0; i < siteVcf.cdr.hdr->nhrec; ++i) {
      bcf_hrec_t *hrec = siteVcf.cdr.hdr->hrec[i];
      //kstring_t s = {0,0,0};
      if ( hrec->type == BCF_HL_FLT ) {
	//hrec_t new_hrec = bcf_hrec_dup(hrec);
	bcf_hdr_add_hrec(odw.hdr,bcf_hrec_dup(hrec));
      }
    } 
  }
      

  // remove INFO tags to be removed
  for(int32_t i=0; i < (int32_t)info2remove.size(); ++i) {
    bcf_hdr_remove(odw.hdr, BCF_HL_INFO, info2remove[i].c_str());
  }
  
  for(int32_t i=0; i < (int32_t)info2update.size(); ++i) {
    bcf_hrec_t* ro = bcf_hdr_get_hrec(odw.hdr, BCF_HL_INFO, "ID", info2update[i].c_str(), NULL);
    bcf_hrec_t* rn = bcf_hdr_get_hrec(siteVcf.cdr.hdr, BCF_HL_INFO, "ID", info2update[i].c_str(), NULL);    
    if ( ro == NULL ) { // then add the INFO field from the site VCF
      if ( rn == NULL )
	error("[E:%s:%d %s] Cannot find INFO field %s from %s", __FILE__, __LINE__, __PRETTY_FUNCTION__, info2update[i].c_str(), siteVcf.bcf_file_name.c_str());
      //notice("foo %s", info2update[i].c_str());
      if ( !bcf_hdr_add_hrec(odw.hdr, bcf_hrec_dup(rn)) )
	error("[E:%s:%d %s] Cannot add INFO field %s from %s to %s", __FILE__, __LINE__, __PRETTY_FUNCTION__, info2update[i].c_str(), siteVcf.bcf_file_name.c_str(), genoVcf.bcf_file_name.c_str());	
    }
    else if ( !same_hrecs(odw.hdr, ro, siteVcf.cdr.hdr, rn) )
      error("[E:%s:%d %s] Incompatible INFO field %s. Add --rm-info %s to avoid errors and completely replace the field with new format", __FILE__, __LINE__, __PRETTY_FUNCTION__, info2update[i].c_str(), info2update[i].c_str());
  }

  bcf_hdr_sync(odw.hdr);
  odw.write_hdr();  

  std::map<std::string,std::string> rsmap;
  if ( !dbsnpVcf.bcf_file_name.empty() ) {
    while( dbsnpVcf.read() ) {
      bcf1_t* dv = dbsnpVcf.cursor();      
      bcf_unpack(dv, BCF_UN_STR);      
      rsmap[dbsnpVcf.get_var_ID()] = dv->d.id;
    }
  }  

  bool sread = false;
  void* values = NULL;
  int32_t nvalues = 0;
  
  while( ( sread = siteVcf.read() ) ) {
    bcf1_t* sv = siteVcf.cursor();

    if ( ( !dbsnpVcf.target_loci.empty() ) && ( !dbsnpVcf.target_loci.overlaps(bcf_get_chrom(siteVcf.cdr.hdr, sv), sv->pos+1, sv->pos+1) ) ) {
      notice("bar %d", sv->pos);      
      continue;
    }

    bool gread = genoVcf.read();
    /*
    while ( gread && ( !genoVcf.target_loci.empty() ) && ( !genoVcf.target_loci.overlaps(bcf_get_chrom(genoVcf.cdr.hdr, sv), sv->pos, sv->pos) ) ) {
      notice("foo %d", sv->pos);
      //continue;
      gread = genoVcf.read();      
    }
    */

    if ( ( !sread ) || (!gread) )
      error("[E:%s:%d %s] Two input VCF files do not share exact site list", __FILE__, __LINE__, __PRETTY_FUNCTION__);
        
    bcf1_t* gv = bcf_dup(genoVcf.cursor());

    bcf_unpack(gv, BCF_UN_SHR);
    bcf_unpack(sv, BCF_UN_SHR);

    if ( ( gv->rid != sv->rid ) || ( gv->pos != sv->pos ) || ( gv->rlen != sv->rlen ) )
      error("[E:%s:%d %s] Two input VCF files do not share exact site list (%d:%d:%d) vs (%d:%d:%d)", __FILE__, __LINE__, __PRETTY_FUNCTION__, gv->rid, gv->pos, gv->rlen, sv->rid, sv->pos, sv->rlen );

    if ( replaceQual )
      gv->qual = sv->qual;

    if ( replaceFilter ) {
      //bcf_update_filter(odw.hdr, gv, sv->d.flt, sv->d.n_flt);
      // Remove all the filters
      for(int32_t i=gv->d.n_flt-1; i >= 0; --i) {
	bcf_remove_filter(odw.hdr, gv, gv->d.flt[i], 0);
      }

      // Add new filters
      for(int32_t i=0; i < sv->d.n_flt; ++i) {
	const char* flt_name = bcf_hdr_int2id(siteVcf.cdr.hdr, BCF_DT_ID, sv->d.flt[i]);
	int32_t new_flt_id = bcf_hdr_id2int(odw.hdr, BCF_DT_ID, flt_name);
	bcf_add_filter(odw.hdr, gv, new_flt_id);
      }      
    }

    if ( !info2remove.empty() ) {
      for(int32_t i=0; i < (int32_t)info2remove.size(); ++i) {
	bcf_update_info(genoVcf.cdr.hdr, gv, info2remove[i].c_str(), NULL, 0, BCF_HT_INT);
      }
    }

    if ( !info2update.empty() ) {
      for(int32_t i=0; i < (int32_t)info2update.size(); ++i) {
	bcf_info_t* info = bcf_get_info(siteVcf.cdr.hdr, sv, info2update[i].c_str());
	if ( info == NULL ) {
	  error("Cannot find %s from site VCF", info2update[i].c_str());
	  bcf_update_info(odw.hdr, gv, info2update[i].c_str(), NULL, 0, BCF_HT_INT);
	}
	else {
	  int32_t htype = 0, ret = 0;
	  switch(info->type) {
	  case BCF_BT_NULL:
	    htype = BCF_HT_FLAG;
	    break;
	  case BCF_BT_INT8: case BCF_BT_INT16: case BCF_BT_INT32:
	    htype = BCF_HT_INT;
	    break;
	  case BCF_BT_FLOAT:
	    htype = BCF_HT_REAL;
	    break;
	  case BCF_BT_CHAR:
	    htype = BCF_HT_STR;
	    break;
	  default:	    
	    error("[E:%s:%d %s] Cannot recognize info->type %d for field %s", __FILE__, __LINE__, __PRETTY_FUNCTION__, info->type, info2update[i].c_str());
	  }
	  if ( (ret = bcf_get_info_values(siteVcf.cdr.hdr, sv, info2update[i].c_str(), &values, &nvalues, htype)) < 0 )
	    error("[E:%s:%d %s] %s: ret=%d, nvalues = %d, type = %d, len = %d, (i,f)=(%d,%f)", __FILE__, __LINE__, __PRETTY_FUNCTION__, info2update[i].c_str(), ret, nvalues, info->type, info->len, info->v1.i, info->v1.f);
	  if ( bcf_update_info(odw.hdr, gv, info2update[i].c_str(), values, nvalues, htype) < 0 )
	    error("[E:%s:%d %s] Failed to update INFO field %s", info2update[i].c_str());
	}
      }      
    }

    /*
    // update variant IDs from dbSNP
    if ( !dbsnpVcf.bcf_file_name.empty() ) {
      const char* chrom = bcf_get_chrom(odw.hdr,gv);
      if ( dbsnpVcf.jump_to(chrom, gv->pos) ) {
	while( dbsnpVcf.read() ) {
	  bcf1_t* dv = dbsnpVcf.cursor();
	  bcf_unpack(dv, BCF_UN_STR);
	  if ( strcmp(chrom,bcf_get_chrom(dbsnpVcf.cdr.hdr, dv)) == 0 ) {
	    if ( dv->pos < gv->pos ) continue;
	    else if ( dv->pos > gv->pos ) break;
	    else if ( dv->rlen != gv->rlen ) continue;
	    else if ( dv->n_allele != gv->n_allele ) continue;
	    else {
	      int32_t i;
	      for(i=0; i < dv->n_allele; ++i) {
		if ( strcmp(dv->d.allele[i], gv->d.allele[i]) != 0 ) break;
	      }
	      if ( i == dv->n_allele ) {
		bcf_update_id(odw.hdr, gv, dv->d.id);
	      }
	      else continue;
	    }
	  }
	  else break;
	}
      }
    }
    */
    if ( !rsmap.empty() ) {
      std::map<std::string,std::string>::iterator it = rsmap.find(siteVcf.get_var_ID());
      if ( it != rsmap.end() )
	bcf_update_id(odw.hdr, gv, it->second.c_str());
    }
    else if ( replaceID ) {
      bcf_update_id(odw.hdr, gv, sv->d.id);      
    }

    //notice("Writing %d:%d", gv->rid, gv->pos);
    odw.write(gv);
    bcf_clear(gv);
    bcf_destroy(gv);
  }
  odw.close();
  //siteVcf.close();
  //genoVcf.close();

  //if ( !dbsnpVcf.bcf_file_name.empty() )
  //dbsnpVcf.close();

  notice("Analysis Finished");

  return 0;
}

