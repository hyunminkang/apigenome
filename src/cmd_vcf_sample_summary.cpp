#include "bcf_filter_arg.h"
#include "cramore.h"
#include "bcf_ordered_reader.h"

int32_t cmdVcfSampleSummary(int32_t argc, char** argv) {
  std::string inVcf;
  std::string out;
  std::string reg;
  std::vector<std::string> sumFields;
  int32_t minDistBp = 0;
  int32_t verbose = 1000;
  bool countVariants = false;  // count variants by variant type, allele type, and genotypes

  bcf_vfilter_arg vfilt;
  bcf_gfilter_arg gfilt;
  
  std::vector<int32_t> acThres;
  std::vector<double> afThres;
  bool posOnly = false;

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input Sites", NULL)
    LONG_STRING_PARAM("in-vcf",&inVcf, "Input VCF/BCF file")
    LONG_STRING_PARAM("region",&reg,"Genomic region to focus on")
    LONG_PARAM("pos-only",&posOnly,"Consider POS field only (not REF field) when determining --region option. This will exclude >1bp variants outside of the region included")

    LONG_PARAM_GROUP("Analysis Options", NULL)
    LONG_MULTI_STRING_PARAM("sum-field",&sumFields, "Field values to calculate the sums")
    LONG_PARAM("count-variants",&countVariants, "Flag to turn on counting variants")    

    LONG_PARAM_GROUP("Variant Filtering Options", NULL)    
    LONG_INT_PARAM("min-dist-bp",&minDistBp, "Minimum distance from the previous variant in base-position")    
    LONG_MULTI_STRING_PARAM("apply-filter",&vfilt.required_filters, "Require at least one of the listed FILTER strings")
    LONG_STRING_PARAM("include-expr",&vfilt.include_expr, "Include sites for which expression is true")
    LONG_STRING_PARAM("exclude-expr",&vfilt.exclude_expr, "Exclude sites for which expression is true")    

    LONG_PARAM_GROUP("Genotype Filtering Options", NULL)    
    LONG_MULTI_INT_PARAM("ac",&acThres,"Allele count threshold to count rare/common variants")
    LONG_MULTI_DOUBLE_PARAM("af",&afThres,"Allele frequency threshold to count rare/common variants")
    LONG_INT_PARAM("minDP",&gfilt.minDP,"Minimum depth threshold for counting genotypes")
    LONG_INT_PARAM("minGQ",&gfilt.minGQ,"Minimum depth threshold for counting genotypes") 

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &out, "Output VCF file name")
    LONG_INT_PARAM("verbose",&verbose,"Frequency of verbose output (1/n)")    
  END_LONG_PARAMS();
  
  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();
  
  // sanity check of input arguments
  if ( inVcf.empty() || out.empty() ) {
    error("[E:%s:%d %s] --in-vcf, --out are required parameters",__FILE__,__LINE__,__FUNCTION__);
  }



  std::vector<GenomeInterval> intervals;
  if ( !reg.empty() ) {
    parse_intervals(intervals, "", reg);
  }
  BCFOrderedReader odr(inVcf, intervals);
  bcf1_t* iv = bcf_init();

  // handle filter string
  std::string filter_str;
  int32_t filter_logic = 0;
  if ( vfilt.include_expr.empty() ) {
    if ( vfilt.exclude_expr.empty() ) {
      // do nothing
    }
    else {
      filter_str = vfilt.exclude_expr;
      filter_logic |= FLT_EXCLUDE;
    }
  }
  else {
    if ( vfilt.exclude_expr.empty() ) {
      filter_str = vfilt.include_expr;
      filter_logic |= FLT_INCLUDE;      
    }
    else {
      error("[E:%s:%d %s] Cannot use both --include-expr and --exclude-expr options",__FILE__,__LINE__,__FUNCTION__);
    }    
  }

  filter_t* filt = NULL;
  if ( filter_logic != 0 )
    filter_init(odr.hdr, filter_str.c_str());

  // handle --apply-filtrs
  std::vector<int32_t> req_flt_ids;
  if ( !vfilt.required_filters.empty() ) {
    for(int32_t i=0; i < (int32_t)vfilt.required_filters.size(); ++i) {
      req_flt_ids.push_back(bcf_hdr_id2int(odr.hdr, BCF_DT_ID, vfilt.required_filters[i].c_str()));
    }
  }  

  notice("Started Reading site information from VCF file");

  std::map< std::string, std::vector<int64_t> > mapFieldSums;
  std::map< std::string, std::vector<int64_t> > mapFieldVars;  
  int32_t nVariant = 0;
  int32_t nsamples = bcf_hdr_nsamples(odr.hdr);

  for(int32_t i=0; i < (int32_t)sumFields.size(); ++i)
    mapFieldSums[sumFields[i]].resize(nsamples, (int64_t)0);

  std::vector<std::string> varFields;

  varFields.push_back("ALL.SNP");
  varFields.push_back("ALL.OTH");
  varFields.push_back("NREF.SNP");
  varFields.push_back("NREF.OTH");  
  varFields.push_back("REF.SNP");
  varFields.push_back("REF.OTH");  
  varFields.push_back("HET.SNP");
  varFields.push_back("HET.OTH");
  varFields.push_back("ALT.SNP");
  varFields.push_back("ALT.OTH");
  varFields.push_back("MISS.SNP");
  varFields.push_back("MISS.OTH");    

  std::sort(acThres.begin(), acThres.end());
  for(int32_t i=0; i < (int32_t)acThres.size(); ++i) {
    char buf[255];
    sprintf(buf, "AC_%d_%d.SNP", i == 0 ? 1 : acThres[i-1]+1, acThres[i]);
    varFields.push_back(buf);
    sprintf(buf, "AC_%d_%d.OTH", i == 0 ? 1 : acThres[i-1]+1, acThres[i]);
    varFields.push_back(buf);    
  }

  std::sort(afThres.begin(), afThres.end());
  for(int32_t i=0; i < (int32_t)afThres.size(); ++i) {
    char buf[255];
    sprintf(buf, "AF_%f_%f.SNP", i == 0 ? 0 : afThres[i-1], afThres[i]);
    varFields.push_back(buf);
    sprintf(buf, "AF_%f_%f.OTH", i == 0 ? 0 : afThres[i-1], afThres[i]);
    varFields.push_back(buf);    
  }
  
  for(int32_t i=0; i < (int32_t)varFields.size(); ++i) {
    mapFieldVars[varFields[i]].resize(nsamples, (int64_t)0);    
  }

  std::vector<int32_t> varMasks(varFields.size());


  int32_t* p_gt = NULL;
  int32_t n_gt = 0;
  
  int32_t* p_fld = NULL;
  int32_t n_fld = 0;
  int32_t prev_rid = -1, prev_pos = -1;
  int32_t nskip = 0;

  int32_t an = 0, ac_alloc = 0, non_ref_ac = 0;
  int32_t* ac = NULL;
    
  for(int32_t k=0; odr.read(iv); ++k) {  // read marker
    if ( k % verbose == 0 )
      notice("Processing %d markers at %s:%d. Skipped %d markers within %d-bp from the previous marker", k, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1, nskip, minDistBp);

    // if minimum distance is specified, skip the variant
    
    if ( ( prev_rid == iv->rid ) && ( iv->pos - prev_pos < minDistBp ) ) {
      ++nskip;
      continue;
    }

    if ( ( !reg.empty() ) && posOnly && ( ( intervals[0].start1 > iv->pos+1 ) || ( intervals[0].end1 < iv->pos+1 ) ) ) {
      notice("With --pos-only option, skipping variant at %s:%d", bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1);      
      ++nskip;
      continue;
    }
    
    bcf_unpack(iv, BCF_UN_FLT);

    // check --apply-filters
    bool has_filter = req_flt_ids.empty() ? true : false;
    if ( ! has_filter ) {
      //notice("%d %d", iv->d.n_flt, (int32_t)req_flt_ids.size());
      for(int32_t i=0; i < iv->d.n_flt; ++i) {
	for(int32_t j=0; j < (int32_t)req_flt_ids.size(); ++j) {
	  if ( req_flt_ids[j] == iv->d.flt[i] )
	    has_filter = true;
	}
      }
    }

    //if ( k % 1000 == 999 ) abort();

    if ( ! has_filter ) { ++nskip; continue; }

    // check filter logic
    if ( filt != NULL ) {
      int32_t ret = filter_test(filt, iv, NULL);
      if ( filter_logic == FLT_INCLUDE ) { if ( !ret)  has_filter = false; }
      else if ( ret ) { has_filter = false; }
    }

    if ( ! has_filter ) { ++nskip; continue; }    

    ++nVariant;

    if ( countVariants ) {
      // calculate AC
      if ( acThres.size() + afThres.size() > 0 ) {
	hts_expand(int, iv->n_allele, ac_alloc, ac);
	an = 0;
	non_ref_ac = 0;
	bcf_calc_ac(odr.hdr, iv, ac, BCF_UN_INFO|BCF_UN_FMT); // get original AC and AN values from INFO field if available, otherwise calculate
	for (int32_t i=1; i<iv->n_allele; i++)
	  non_ref_ac += ac[i];
	for (int32_t i=0; i<iv->n_allele; i++)
	  an += ac[i];      
      }
      
      // determine variant type : flags are   MISSING HOMALT HET HOMREF 
      std::fill(varMasks.begin(), varMasks.end(), 0);
      if ( bcf_is_snp(iv) ) {
	varMasks[0] = MASK_GT_ALL;
	varMasks[2] = MASK_GT_NONREF;
	varMasks[4] = MASK_GT_HOMREF;
	varMasks[6] = MASK_GT_HET;
	varMasks[8] = MASK_GT_HOMALT;
	varMasks[10] = MASK_GT_MISS;
      } // for non-ref genotypes
      else {
	varMasks[1] = MASK_GT_ALL;
	varMasks[3] = MASK_GT_NONREF;
	varMasks[5] = MASK_GT_HOMREF;
	varMasks[7] = MASK_GT_HET;
	varMasks[9] = MASK_GT_HOMALT;
	varMasks[11] = MASK_GT_MISS;
      } // for non-ref genotypes
      for(int32_t i=0, j=12; i < (int32_t)acThres.size(); ++i, j += 2) {
	if ( ( non_ref_ac > (i == 0 ? 0 : acThres[i-1]) ) && ( non_ref_ac <= acThres[i] ) ) {
	  if ( bcf_is_snp(iv) ) {
	    varMasks[j] = MASK_GT_NONREF;
	  }
	  else {
	    varMasks[j+1] = MASK_GT_NONREF;	    
	  }
	}
      }
      for(int32_t i=0, j=12 + 2*acThres.size(); i < (int32_t)afThres.size(); ++i, j += 2) {
	if ( ( non_ref_ac > (i == 0 ? 0 : afThres[i-1]*an) ) && ( non_ref_ac <= afThres[i]*an ) ) {
	  if ( bcf_is_snp(iv) ) {
	    varMasks[j] = MASK_GT_NONREF;
	  }
	  else {
	    varMasks[j+1] = MASK_GT_NONREF;	    
	  }	  
	}
      }
      
      // extract genotype and apply genotype level filter
      if ( bcf_get_genotypes(odr.hdr, iv, &p_gt, &n_gt) < 0 ) {
	error("[E:%s:%d %s] Cannot find the field GT from the VCF file at position %s:%d",__FILE__,__LINE__,__FUNCTION__, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1);
      }
      for(int32_t i=0; i < nsamples; ++i) {
	int32_t g1 = p_gt[2*i];
	int32_t g2 = p_gt[2*i+1];
	int32_t geno;
	if ( bcf_gt_is_missing(g1) || bcf_gt_is_missing(g2) ) {
	  geno = 0;
	}
	else {
	  geno = ((bcf_gt_allele(g1) > 0) ? 1 : 0) + ((bcf_gt_allele(g2) > 0) ? 1 : 0) + 1;
	  //if ( i == 0 )
	  //notice("g1 = %d, g2 = %d, geno = %d", g1,g2,geno);
	}
	p_gt[i] = geno;
      }

      if ( gfilt.minDP > 0 ) {
	if ( bcf_get_format_int32(odr.hdr, iv, "DP", &p_fld, &n_fld) < 0 ) {
	  error("[E:%s:%d %s] Cannot find the field DP from the VCF file at position %s:%d",__FILE__,__LINE__,__FUNCTION__, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1);

	}
	for(int32_t i=0; i < nsamples; ++i) {
	  if ( p_fld[i] < gfilt.minDP ) p_gt[i] = 0;
	}
      }

      if ( gfilt.minGQ > 0 ) {
	if ( bcf_get_format_int32(odr.hdr, iv, "GQ", &p_fld, &n_fld) < 0 ) {
	  error("[E:%s:%d %s] Cannot find the field GQ from the VCF file at position %s:%d",__FILE__,__LINE__,__FUNCTION__, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1);

	}
	for(int32_t i=0; i < nsamples; ++i) {
	  if ( p_fld[i] < gfilt.minGQ ) p_gt[i] = 0;
	}
      }

      // update the maps
      for(int32_t i=0; i < (int32_t)varFields.size(); ++i) {
	std::vector<int64_t>& v = mapFieldVars[varFields[i]];
	for(int32_t j=0; j < nsamples; ++j) {
	  if ( varMasks[i] & ( 0x01 << p_gt[j] ) ) {
	    //if ( j == 0 ) notice("VarMask[i] = %d, p_gt[j] = %d", varMasks[i], p_gt[j]);
	    ++v[j];
	  }
	}
      }
      //if ( rand() % 100 == 0 ) abort();
             
    }
    
    // perform sumField tasks
    for(int32_t i=0; i < (int32_t)sumFields.size(); ++i) {
      if ( bcf_get_format_int32(odr.hdr, iv, sumFields[i].c_str(), &p_fld, &n_fld) < 0 ) {
	error("[E:%s:%d %s] Cannot find the field %s from the VCF file at position %s:%d",__FILE__,__LINE__,__FUNCTION__, sumFields[i].c_str(), bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1);
      }
      if ( nsamples != n_fld )
	error("[E:%s:%d %s] Field %s has multiple elements",__FILE__,__LINE__,__FUNCTION__,sumFields[i].c_str());
      std::vector<int64_t>& v = mapFieldSums[sumFields[i]];
      if ( (int32_t)v.size() != nsamples )
	error("[E:%s:%d %s] mapFieldSums object does not have %s as key",__FILE__,__LINE__,__FUNCTION__,sumFields[i].c_str());

      for(int32_t j=0; j < nsamples; ++j) {
	//if ( p_fld[j] != bcf_int32_missing ) {
	v[j] += p_fld[j];
	//}
      }
    }

    prev_rid = iv->rid;
    prev_pos = iv->pos;
  }

  htsFile* wf = hts_open(out.c_str(), "w");
  
  hprintf(wf, "ID\tN.VAR");
  for(int32_t i=0; i < (int32_t)varFields.size(); ++i) {
    hprintf(wf, "\t%s",varFields[i].c_str());    
  }
  for(int32_t i=0; i < (int32_t)sumFields.size(); ++i) {
    hprintf(wf, "\tSUM.%s",sumFields[i].c_str());
  }
  hprintf(wf,"\n");
  
  for(int32_t i=0; i < nsamples; ++i) {
    hprintf(wf, "%s", odr.hdr->id[BCF_DT_SAMPLE][i]);
    hprintf(wf, "\t%d", nVariant);
    for(int32_t j=0; j < (int32_t)varFields.size(); ++j) {
      hprintf(wf, "\t%lld", mapFieldVars[varFields[j]][i]);
    }    
    for(int32_t j=0; j < (int32_t)sumFields.size(); ++j) {
      hprintf(wf, "\t%lld", mapFieldSums[sumFields[j]][i]);
    }
    hprintf(wf, "\n");
  }
  hts_close(wf);
  odr.close();

  return 0;
}

