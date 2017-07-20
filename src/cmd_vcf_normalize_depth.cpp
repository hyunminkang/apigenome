#include "cramore.h"
#include "bcf_filtered_reader.h"
#include "bcf_ordered_reader.h"
#include "sex_ploidy_map.h"
#include "fastaGC.h"

int32_t cmdVcfNormalizeDepth(int32_t argc, char** argv) {
  BCFFilteredReader bfr;
  std::string gcFile;
  std::string outPrefix;
  std::string knownVcf;
  bool xyFlag = false;
  bool includeIndels = false;
  std::string xLabel("X");
  std::string yLabel("Y");
  std::string mtLabel("MT");
  int32_t xStart = 2699520;
  int32_t xStop = 154931044;
  int32_t minDist = 100;
  int32_t capDepth = 100;
  int32_t seed = 0;

  bfr.vfilt.maxAlleles = 2;
  bfr.verbose = 1000000;  

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Running modes",NULL)
    LONG_PARAM("xy",&xyFlag,"Compute normalized depths for X/Y chromosomes")

    LONG_PARAM_GROUP("Common input arguments", NULL)
    LONG_STRING_PARAM("vcf",&bfr.bcf_file_name, "Input VCF/BCF filename")
    LONG_STRING_PARAM("known",&knownVcf, "Known sites to consider (optional)")    
    LONG_STRING_PARAM("gc",&gcFile, "GC content file name")
    LONG_STRING_PARAM("out",&outPrefix,"Output prefix")
    LONG_INT_PARAM("seed",&seed,"Random seed")    

    LONG_PARAM_GROUP("Arguments for --xy mode", NULL)
    LONG_PARAM("include-indels",&includeIndels,"Include indels")
    LONG_INT_PARAM("min-dist",&minDist,"Minimum distance between SNPs (within a cluster a SNP is randomly chosen)")
    LONG_INT_PARAM("cap-depth",&capDepth,"Maximum value of depth to cap to")
    
    LONG_STRING_PARAM("xLabel",&xLabel, "Label for X chromosome : X for GRCh37, chrX for GRCh38")
    LONG_STRING_PARAM("yLabel",&yLabel, "Label for Y chromosome : Y for GRCh37, chrY for GRCh38")
    LONG_STRING_PARAM("mtLabel",&mtLabel, "Label for MT chromosome : MT for GRCh37, chrM for GRCh38")
    LONG_INT_PARAM("xStart",&xStart,"Start base position for non-PAR X chromosome : 2699520 for GRCh37, and 2781479 for GRCh38")
    LONG_INT_PARAM("xStop",&xStop,"Stop base position for non-PAR X chromsoomes : 1549311044 for GRCh37, and 155701383 for GRCh38")
  END_LONG_PARAMS();
  
  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  if ( !includeIndels ) bfr.vfilt.snpOnly = true;

  notice("Reading in BCFs...");
  bfr.init_params();

  std::vector< std::set<int32_t> > knownPos;
  if ( !knownVcf.empty() ) {
    notice("Reading in Known Position BCF/VCFs...");
    std::vector<GenomeInterval> intervals;
    BCFOrderedReader odr(knownVcf.c_str(), intervals);
    bcf1_t* v = bcf_init1();
    while( odr.read(v) ) {
      int32_t tid = bcf_hdr_name2id(bfr.cdr.hdr, bcf_hdr_id2name(odr.hdr,v->rid));
      if ( tid >= (int32_t) knownPos.size() ) knownPos.resize(tid+1);
      if ( tid < 0 )
	error("Cannot find chromosome %s at %d", bcf_hdr_id2name(odr.hdr,v->rid), v->rid);
      knownPos[tid].insert(v->pos);
    }
    bcf_destroy(v);
  }
  
  notice("Loading GC profile.. %s", gcFile.c_str());
  fastaGC fGC;
  fGC.on_memory = true;
  if ( !fGC.openGC(gcFile.c_str()) )
    error("Failed loading fastaGC file %s",gcFile.c_str());

  notice("Finished loading GC profile.. %s", gcFile.c_str());  

  // infer autosomes and X Y chromosomes
  std::vector<int32_t> auto_tids;
  int32_t x_tid = -1;
  int32_t y_tid = -1;
  //int32_t mt_tid = -1;
  for(int32_t i=0; i < (int32_t)fGC.seqnames.size(); ++i) {
    if ( fGC.seqnames[i].compare(xLabel) == 0 )       x_tid = i;
    else if ( fGC.seqnames[i].compare(yLabel) == 0 )  y_tid = i;
    else if ( fGC.seqnames[i].compare(mtLabel) == 0 ) {} // mt_tid = i;
    else {
      std::string core = ( fGC.seqnames[i].size() > 3 && fGC.seqnames[i].compare(0,3,"chr") == 0 ) ? fGC.seqnames[i].substr(3) : fGC.seqnames[i];
      //if ( i < 30 ) notice("%s | %s",fGC.seqnames[i].c_str(), core.c_str());      
      int32_t j=0;
      for(j=0; j < (int32_t) core.size(); ++j) {
	if ( ! isdigit(core[j]) ) break;
      }
      if ( j == (int32_t)core.size() ) auto_tids.push_back(i);
    }
  }

  notice("Identified %u autosomes", auto_tids.size());

  if ( xyFlag ) {
    if ( gcFile.empty() || outPrefix.empty() )
      error("--gc and --out required with --xy");
    notice("Calculating normalized XY depths");
    if ( ( x_tid < 0 ) || ( y_tid < 0 ) )
      error("Cannot find X and/or Y chromosome labels from GC content file %s", gcFile.c_str());

    // translate tid from BCF to the AUTO/X/Y/OTHER categories

    std::map<int32_t,int32_t> bcf_tid2type;
    for(int32_t i=0; i < (int32_t)auto_tids.size(); ++i) {
      int32_t tid = bcf_hdr_name2id(bfr.cdr.hdr, fGC.seqnames[auto_tids[i]].c_str());
      if ( tid < 0 ) error("Cannot find chromosome %s from BCF header",fGC.seqnames[auto_tids[i]].c_str());
      bcf_tid2type[tid] = PLOIDY_TYPE_AUTOSOME;
      //notice("Autosome: %s, %d",fGC.seqnames[auto_tids[i]].c_str(), auto_tids[i]);      
    }
    int32_t tid = bcf_hdr_name2id(bfr.cdr.hdr, fGC.seqnames[x_tid].c_str());
    if ( tid < 0 ) error("Cannot find chromosome %s from BCF header",fGC.seqnames[x_tid].c_str());
    bcf_tid2type[tid] = PLOIDY_TYPE_X;

    tid = bcf_hdr_name2id(bfr.cdr.hdr, fGC.seqnames[y_tid].c_str());
    if ( tid < 0 ) error("Cannot find chromosome %s from BCF header",fGC.seqnames[y_tid].c_str());    
    bcf_tid2type[tid] = PLOIDY_TYPE_Y;    

    int32_t cur_rid = -1;
    int32_t cur_type = -1;
    int32_t cur_beg = -1;
    int32_t cur_end = -1;
    int32_t cur_size = -1;

    // store histogram of (GC, DEPTH) as uint32_t;
    typedef std::map<uint32_t,int64_t> gc_dp_hist_t;
    typedef std::map<uint32_t,int64_t>::iterator gc_dp_hist_itr_t;
    
    std::vector<gc_dp_hist_t> hists;
    //gc_dp_hist_t cur_hist;
    uint32_t cur_clust = 0xffffffff;
    
    if ( hists.size() < PLOIDY_TYPE_AUTOSOME + 1 )
      hists.resize(PLOIDY_TYPE_AUTOSOME+1);
    if ( hists.size() < PLOIDY_TYPE_X + 1 )
      hists.resize(PLOIDY_TYPE_X+1);
    if ( hists.size() < PLOIDY_TYPE_Y + 1 )
      hists.resize(PLOIDY_TYPE_Y+1);

    int32_t* dps = NULL;
    int32_t n_dps = 0;
    int64_t nsnps = 0;
    int32_t n_types[4] = {0,0,0,0};
    
    while( bfr.read() ) {
      bcf1_t* v = bfr.cursor();

      if ( (!includeIndels) && ( ! bcf_is_snp(v) ) ) continue;

      if ( ( ! knownPos.empty() ) && ( ( (int32_t)knownPos.size() <= v->rid ) || ( knownPos[v->rid].find(v->pos) == knownPos[v->rid].end() ) ) ) continue;

      // determine the ploidy
      if ( cur_rid != v->rid ) {
	std::map<int32_t,int32_t>::iterator it = bcf_tid2type.find(v->rid);
	if ( it == bcf_tid2type.end() ) continue;

	cur_rid = v->rid;
	cur_type = it->second;
	cur_beg = -1;
	cur_end = -1;
	cur_size = -1;
	cur_clust = 0xffffffff;
      }

      if ( ( cur_type == PLOIDY_TYPE_X ) && ( ( v->pos < xStart ) || ( v->pos >= xStop ) ) ) continue;      

      uint16_t gc = fGC.getGC(bcf_hdr_id2name(bfr.cdr.hdr,v->rid), v->pos+1);
      bcf_unpack(v, BCF_UN_ALL);

      if ( bcf_get_format_int32(bfr.cdr.hdr, v, "N", &dps, &n_dps) < 0 )
	error("Cannot find FORMAT field N at %s:%d", bcf_hdr_id2name(bfr.cdr.hdr,v->rid), v->pos+1);

      if ( n_dps != 1 )
	error("More than one FORMAT fields N observed at %s:%d", bcf_hdr_id2name(bfr.cdr.hdr,v->rid), v->pos+1);

      if ( dps[0] > capDepth ) dps[0] = capDepth;

      uint32_t gcdp = ( ( gc << 16 ) | ((uint32_t)dps[0]) );

      if ( ( cur_beg < 0 ) || ( v->pos - cur_end >= minDist ) ) {
	cur_beg = v->pos;
	cur_end = v->pos;
	cur_size = 1;
	cur_clust = gcdp;
	
	// need to insert new
	++(hists[cur_type][gcdp]);
	++nsnps;
	++n_types[cur_type];
      }
      else { // clustered
	cur_end = v->pos;
	++cur_size;
	double p = (rand()+0.5) / (RAND_MAX + 1.0);
	if ( p * cur_size < 1 ) { // need to replace
	  --(hists[cur_type][cur_clust]);
	  ++(hists[cur_type][gcdp]);	  
	  cur_clust = gcdp;
	}
      }
    }

    //bfr.close();

    notice("Finished reading %lld variants from BCF/VCF files to get GC-adjusted histograms", nsnps);

    // Normalize the depths carefully
    // Calculate autosomes depth first
    gc_dp_hist_t& h = hists[PLOIDY_TYPE_AUTOSOME];
    std::map<uint16_t,double> gc_fracs;

    double gc_sum = 0;
    double dp_sum = 0;
    // it->first represents GC/DP it->second represents # sites
    for(gc_dp_hist_itr_t it = h.begin(); it != h.end(); ++it) {
      uint16_t gc = (uint16_t)(( it->first >> 16 ) & 0x0ffff);
      uint16_t dp = (uint16_t)(it->first & 0x0ffff);
      if ( gc == (uint16_t)0xffff ) continue;

      // gc_fracs will contain total number of variants for each depth
      gc_fracs[gc] += it->second;

      if ( dp > 100 ) { error("dp = %u, gc = %u", (uint32_t)dp, (uint32_t)gc); }
      
      gc_sum += (double) it->second;
      dp_sum += ((double) dp * (double) it->second );
    }

    double auto_depth = dp_sum / gc_sum;
    for(std::map<uint16_t,double>::iterator it = gc_fracs.begin();
	it != gc_fracs.end(); ++it) {
      it->second /= gc_sum;
    }

    gc_sum = 0;
    h = hists[PLOIDY_TYPE_X];
    
    std::map<uint16_t,double> x_fracs;
    std::map<uint16_t,double> x_dpsum;
    
    for(gc_dp_hist_itr_t it = h.begin(); it != h.end(); ++it) {
      uint16_t gc = (uint16_t)(( it->first >> 16 ) & 0x0ffff);
      uint16_t dp = (uint16_t)(it->first & 0x0ffff);
      if ( gc == (uint16_t)0xffff ) continue;
      if ( dp > 100 ) { error("dp = %u", dp); }

      if ( it->second == 0 ) continue; // error("it = (%u,%u) => (%lg)", (uint32_t)gc, (uint32_t)dp, it->second);

      x_fracs[gc] += (double) it->second;
      x_dpsum[gc] += ( (double) dp * (double) it->second );

      if ( x_dpsum[gc] < 0 )
	error("%u %u %lld %lg %lg", gc, dp, it->second, x_fracs[gc], x_dpsum[gc]);

      gc_sum += (double)it->second; // total number of variants
    }

    double x_depth = 0;
    for(std::map<uint16_t,double>::iterator it = x_fracs.begin();
	it != x_fracs.end(); ++it) {
      it->second /= gc_sum;
      x_depth += (x_dpsum[it->first] / gc_sum * gc_fracs[it->first] / it->second);
    }

    gc_sum = 0;
    h = hists[PLOIDY_TYPE_Y];
    
    std::map<uint16_t,double> y_fracs;
    std::map<uint16_t,double> y_dpsum;
    
    for(gc_dp_hist_itr_t it = h.begin(); it != h.end(); ++it) {
      uint16_t gc = (uint16_t)(( it->first >> 16 ) & 0x0ffff);
      uint16_t dp = (uint16_t)(it->first & 0x0ffff);
      if ( gc == (uint16_t)0xffff ) continue;
      if ( dp > 100 ) { error("dp = %u", dp); }

      if ( it->second == 0 ) continue;

      y_fracs[gc] += (double) it->second;
      y_dpsum[gc] += ( (double) dp * (double) it->second );

      gc_sum += (double)it->second;
    }

    double y_depth = 0;
    for(std::map<uint16_t,double>::iterator it = y_fracs.begin();
	it != y_fracs.end(); ++it) {
      it->second /= gc_sum;
      y_depth += ((y_dpsum[it->first] / gc_sum) * gc_fracs[it->first] / it->second);
    }

    notice("AUTO_DEPTH = %.2lf", auto_depth);
    notice("AUTO_COUNT = %d",    n_types[PLOIDY_TYPE_AUTOSOME]);    
    notice("X_DEPTH = %.2lf",    x_depth);
    notice("X_COUNT = %d",       n_types[PLOIDY_TYPE_X]);        
    notice("Y_DEPTH = %.2lf",    y_depth);
    notice("Y_COUNT = %d",       n_types[PLOIDY_TYPE_Y]);            
    notice("X_REL_DEPTH = %.2lf", x_depth / auto_depth);
    notice("Y_REL_DEPTH = %.2lf", y_depth / auto_depth);

    htsFile* ofh = hts_open((outPrefix+".xy").c_str(), "w");
    hprintf(ofh,"CHRX_REL_DEPTH\t%.3lf\n",x_depth/auto_depth);
    hprintf(ofh,"CHRY_REL_DEPTH\t%.3lf\n",y_depth/auto_depth);  
    hprintf(ofh,"AUTO_DEPTH\t%.3lf\n",auto_depth);
    hprintf(ofh,"CHRX_DEPTH\t%.3lf\n",x_depth);
    hprintf(ofh,"CHRY_DEPTH\t%.3lf\n",y_depth);
    hprintf(ofh,"AUTO_COUNT\t%.2lf\n",n_types[PLOIDY_TYPE_AUTOSOME]);
    hprintf(ofh,"CHRX_COUNT\t%.2lf\n",n_types[PLOIDY_TYPE_X]);
    hprintf(ofh,"CHRY_COUNT\t%.2lf\n",n_types[PLOIDY_TYPE_Y]);
    hts_close(ofh);
  }

  return 0;
}
