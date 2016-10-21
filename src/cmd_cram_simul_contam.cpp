#include "cramore.h"

KHASH_MAP_INIT_STR(vdict, bcf_idinfo_t)
typedef khash_t(vdict) vdict_t;

int32_t cmdCramSimulContam(int32_t argc, char** argv) {
  std::string inSam;
  std::string inVcf;
  std::string inContamMap;
  std::string outPrefix;
  std::string tagGroup;
  std::string tagUMI;
  std::vector<std::string> smContams;
  int32_t capBQ = 40;
  int32_t seed = 0;
  int32_t verbose = 10000;

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Options for input SAM/BAM/CRAM", NULL)
    LONG_STRING_PARAM("sam",&inSam, "Input SAM/BAM/CRAM file. Must be sorted by coordinates and indexed")
    LONG_STRING_PARAM("tag-group",&tagGroup, "Tag representing readgroup or cell barcodes, in the case to partition the BAM file into multiple groups")
    LONG_STRING_PARAM("tag-UMI",&tagUMI, "Tag representing UMIs")

    LONG_PARAM_GROUP("Options for input VCF/BCF", NULL)
    LONG_STRING_PARAM("vcf",&inVcf, "Input VCF/BCF file")

    LONG_PARAM_GROUP("Options for specifying distribution of contamination", NULL)
    LONG_STRING_PARAM("contam-map",&inContamMap, "Each line contains [BARCODE] [ID1,FRAC_READS_1] [ID2,FRAC_READ_2] ...")
    LONG_MULTI_STRING_PARAM("sm-contam",&smContams, "List of string representing [SAMPLE_ID,FRAC_READS_1]")
    LONG_INT_PARAM("seed",&seed,"Randomization seed")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out",&outPrefix,"Out prefix")

    LONG_PARAM_GROUP("Analysis Options", NULL)
    LONG_INT_PARAM("cap-BQ", &capBQ, "Maximum base quality (higher BQ will be capped)")
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
  samFile* out = NULL;

  if ( ( in = sam_open(inSam.c_str(), "r") ) == 0 ) {
    error("Cannot open file %s for reading\n",inSam.c_str());    
  }

  if ( ( out = sam_open(outPrefix.c_str(), "w") ) == 0 ) {
    error("Cannot open file %s for writing\n",outPrefix.c_str());    
  }  

  if ( ( header = sam_hdr_read(in) ) == 0 ) {
    error("Cannot open header from %s\n",inSam.c_str());
  }

  if ( sam_hdr_write(out, header) < 0 ) {
    error("Cannot write header to %s\n",outPrefix.c_str());    
  }

  if ( seed == 0 )
    srand(time(NULL));
  else
    srand(seed);

  bam1_t *b = bam_init1();
  //int32_t r;    

  int32_t nsamples = bcf_hdr_nsamples(odr.hdr);
  std::vector< std::vector<uint8_t> > v_gts;
  std::vector<double> v_afs;
  std::vector<int32_t> v_rids;
  std::vector<int32_t> v_poss;
  std::vector<char> v_refs;
  std::vector<char> v_alts;
  
  // identify samples to focus on
  std::map<std::string, int32_t> sm2idx;
  for(int32_t i=0; i < nsamples; ++i) {
    sm2idx[odr.hdr->samples[i]] = i;
  }

  // load contaminaton map
  std::map<std::string, std::vector< std::pair<int32_t, double> > > mapContam;
  
  if ( inContamMap.empty() ) {
    if ( smContams.empty() ) {
      error("Either --contan-map or --sm-contam arguments are required");
    }
    std::vector< std::pair<int32_t, double> >& v = mapContam["."]; // default barcode
    for(int32_t i=0; i < (int32_t)smContams.size(); ++i) {
      uint32_t icomma = smContams[i].find(',');
      if ( icomma == std::string::npos ) {
	error("Cannot recognize --sm-contam %s. Must contain comma",smContams[i].c_str());
      }
      std::string id = smContams[i].substr(0, icomma);
      double alpha = atof(smContams[i].substr(icomma+1).c_str());
      if ( sm2idx.find(id) == sm2idx.end() )
	error("Cannot find sample ID %s from the VCF",id.c_str());
      v.push_back( std::pair<int32_t, double>(sm2idx[id], alpha) );
    }
  }
  else {
    htsFile* hp = hts_open(inContamMap.c_str(), "r");
    if ( hp == NULL )
      error("Cannot open file %s for reading", inContamMap.c_str());

    kstring_t str = {0,0,0};
    int32_t lstr;
    int32_t nfields = 0;
    int32_t* fields = NULL;
    while( ( lstr = hts_getline(hp, KS_SEP_LINE, &str) ) >= 0 ) {
      fields = ksplit(&str, 0, &nfields);
      if ( nfields < 2 )
	error("Cannot parse line %s to extract contamination fraction",str.s);

      std::vector< std::pair<int32_t, double> >& v = mapContam[&str.s[fields[0]]];
      
      for(int32_t i=1; i < nfields; ++i) {
	char* pcomma = strchr(&str.s[fields[i]], ',');
	if ( pcomma == NULL ) {
	  error("Cannot recognize %s. Must contain comma",&str.s[fields[i]]);  
	}
	std::string id = std::string(&str.s[fields[i]], pcomma-&str.s[fields[i]]);
	double alpha = atof(pcomma+1);
	if ( sm2idx.find(id) == sm2idx.end() )
	  error("Cannot find sample ID %s from the VCF",id.c_str());
	v.push_back( std::pair<int32_t, double>(sm2idx[id], alpha) );
      }
    }
  }

  // normalize the contamination fraction
  for(std::map<std::string, std::vector< std::pair<int32_t, double> > >::iterator it = mapContam.begin(); it != mapContam.end(); ++it) {
    std::vector< std::pair<int32_t, double > >& v = it->second;
    double sumAlpha = 0;
    for(int32_t i=0; i < (int32_t)v.size(); ++i) {
      sumAlpha += v[i].second;
    }
    if ( fabs(sumAlpha - 1) > 1e-8 )
      notice("Sum of alphas for barcode %s is %lg, normalizing to be 1..", it->first.c_str(), sumAlpha);

    for(int32_t i=0; i < (int32_t)v.size(); ++i) {
      v[i].second /= sumAlpha;
    }    
  }

  notice("Started reading from the VCF file %s", inVcf.c_str());
  
  // read VCF and store genotypes
  while( odr.read(iv) ) { // read marker
    //bool skip = false;
    bcf_unpack(iv, BCF_UN_ALL);
    
    if ( iv->n_allele > 2 ) continue; // skip multi-allelics
    if ( !bcf_is_snp(iv) ) continue;  // focus only on SNPs
    
    // chrom is iv->rid
    // position is iv->pos
    // ref is iv->d.allele[0]
    // read genotypes
    int32_t rid = iv->rid;
    int32_t pos = iv->pos;
    char ref = iv->d.allele[0][0];
    char alt = iv->d.allele[1][0];

    // store marker information
    v_rids.push_back(rid);
    v_poss.push_back(pos);
    v_refs.push_back(ref);
    v_alts.push_back(alt);
    
    uint32_t* p_gts = (uint32_t*)calloc(nsamples * 2, sizeof(uint32_t));
    int32_t n_gts = 0;
    
    // extract genotypes fpr selected individuals
    if ( bcf_get_genotypes(odr.hdr, iv, &p_gts, &n_gts) < 0 ) {
      error("Cannot extract genotypes at %s:%d %c/%c", bcf_hdr_id2name(odr.hdr, rid), pos, ref, alt);
    }
    if ( n_gts != nsamples * 2 ) {
      error("Cannot extract %d genotypes at %s:%d %c/%c. Extracted only %d", nsamples * 2, bcf_hdr_id2name(odr.hdr, rid), pos, ref, alt, n_gts); 
    }

    v_gts.resize( v_gts.size() + 1 );
    //v_gts.push_back( std::vector<uint8_t>(nv, 0) );
    std::vector<uint8_t>& v = v_gts.back();
    v.resize(nsamples);

    int32_t ac = 0, an = 0;
    for(int32_t i=0; i < nsamples; ++i) {   // bi-allelic encoding of variant
      int32_t g1 = p_gts[2*i];
      int32_t g2 = p_gts[2*i+1];
      uint8_t geno;
      if ( bcf_gt_is_missing(g1) || bcf_gt_is_missing(g2) ) {
	geno = 0;
      }
      else {
	geno = (uint8_t)(((bcf_gt_allele(g1) > 0) ? 1 : 0) + ((bcf_gt_allele(g2) > 0) ? 1 : 0) + 1);
	ac += (geno-1);
	an += 2;
      }
      v[i] = geno;
    }

    if ( ( verbose > 0 ) && ( v_rids.size() % verbose == 0 ) )
      notice("Reading %d markers from the VCF file", (int32_t)v_rids.size());

    v_afs.push_back((double)(ac+5e-7)/(double)(an+1e-6));

    free(p_gts);
  }

  notice("Finished reading %d markers from the VCF file", (int32_t)v_rids.size());

  std::map<int32_t, int32_t> rid_vcf2sam;
  vdict_t *d = (vdict_t*)odr.hdr->dict[BCF_DT_CTG];
  for(khint_t k=kh_begin(d); k < kh_end(d); ++k) {
    if ( !kh_exist(d,k) ) continue;
    int32_t vtid = kh_val(d,k).id;
    int32_t btid = bam_name2id(header, kh_key(d,k));
    if ( btid >= 0 ) 
      rid_vcf2sam[vtid] = btid;
  }

  // ensure that the chromosome orders are increasing
  int32_t prev = -1;
  for(std::map<int32_t,int32_t>::const_iterator it = rid_vcf2sam.begin(); it != rid_vcf2sam.end(); ++it) {
    if ( prev >= it->second )
      error("The references sequences are not ordered consistently between BAM and VCF file");
    prev = it->second;
  }
  
  // reading from BAM files
  int32_t ibeg = -1;
  int32_t iend = -1;
  int32_t ichr = -1;

  char gtag[2] = {0,0};
  //char utag[2] = {0,0};
  if ( tagGroup.empty() ) { // do nothing
  }
  else if ( tagGroup.size() == 2 ) {
    gtag[0] = tagGroup.at(0);
    gtag[1] = tagGroup.at(1);    
  }
  else {
    error("Cannot recognize group tag %s. It is suppose to be a length 2 string",tagGroup.c_str());
  }

  if ( tagUMI.empty() ) { // do nothing
  }
  else if ( tagUMI.size() == 2 ) {
    //utag[0] = tagUMI.at(0);
    //utag[1] = tagUMI.at(1);    
  }
  else {
    error("Cannot recognize UMI tag %s. It is suppose to be a length 2 string",tagUMI.c_str());
  }

  int32_t nReadsKept = 0, nReadsChanged = 0;
  int32_t ret = 0;
  kstring_t readseq = {0,0,0};
  kstring_t readqual = {0,0,0};
  char base, qual;
  int32_t rpos;

  while( ( ret = sam_read1(in, header, b) ) >= 0 ) {
    // check whether any variant overlaps with the reads
    int32_t tid = b->core.tid;
    int32_t beg1 = bam_get_pos1(b);
    int32_t end1 = bam_get_end_pos1(b);

    // advance indices
    while( ( ichr < tid ) || ( ( ichr == tid ) && ( v_poss[ibeg] < beg1 ) ) ) {
      ++ibeg;
      ichr = (rid_vcf2sam.find(v_rids[ibeg]) == rid_vcf2sam.end() ? -1 : rid_vcf2sam[v_rids[ibeg]]);
    }

    if ( ichr == tid ) {
      if ( v_poss[ibeg] < beg1 ) {
	++nReadsKept;
	// do nothing, just print it
	sam_write1(out, header, b);
      }
      else {
	iend = ibeg;
	while( ( v_rids[iend] == v_rids[ibeg] ) && ( v_poss[iend] <= end1 ) ) {
	  ++iend;
	}

	// determine the originating sample first
	uint8_t *bcd = (*gtag) ? (uint8_t*) bam_aux_get(b, gtag) : NULL;
	std::vector< std::pair<int32_t,double> >& v = mapContam["."];
	if ( bcd != NULL ) {
	  if ( *bcd == 'Z' ) {
	    const char* sbcd = bam_aux2Z(bcd);
	    v = mapContam[sbcd];
	    if ( v.size() == 0 ) {
	      notice("Cannot find barcode %s. Skippingg...",sbcd);
	      continue;
	    }
	  }
	}

	if ( v.size() == 0 ) {
	  notice("Cannot find barcode . Skippingg...");
	  continue;	  
	}

	// randomly sample the originting sample
	double r = (rand()+0.5)/(RAND_MAX+1.);
	double ir = v.size()-1;
	for(int32_t i=0; i < (int32_t)v.size()-1; ++i) {
	  if ( r < v[i].second ) {
	    ir = i;
	    break;
	  }
	  r -= v[i].second;
	}


	// modify corresponding bases
	for(int32_t i=ibeg; i < iend; ++i) {
	  uint16_t iref = seq_nt16_table[(int32_t)v_refs[i]];
	  uint16_t ialt = seq_nt16_table[(int32_t)v_alts[i]];	
	  
	  // modify corresponding bases	  
	  bam_get_base_and_qual_and_read_and_qual(b, (uint32_t)v_poss[i], base, qual, rpos, &readseq, &readqual);

	  if ( rpos >= 0 ) {
	    r = (rand()+0.5)/(RAND_MAX+1.);		
	    
	    uint8_t gt = v_gts[i][v[ir].first];
	    bool ref = false;
	    if ( gt == 0 ) {
	      ref = (r < v_afs[i]) ? false : true;
	    }
	    else if ( gt == 1 ) {
	      ref = true;
	    }
	    else if ( gt == 2 ) {
	      ref = ( r < 0.5 ) ? false : true;
	    }
	    else if ( gt == 3 ) {
	      ref = false;
	    }

	    r = (rand()+0.5)/(RAND_MAX+1.);
	    if ( r < phredConv.phred2Err[readqual.s[rpos]-33] ) {
	      ref = ( ref ? false : true );
	    }
	    bam_get_seq(b)[rpos] = (ref ? iref : ialt);

	    notice("Modified: %d %d %d %d\n", ichr, v_poss[i], rpos, bam_get_seq(b)[rpos]);
	  }
	}
	sam_write1(out, header, b);	  
      }
    }
  }

  notice("Finished writing BAM file. %d reads unchanged, %d reads changed", nReadsKept, nReadsChanged);

  odr.close();
  hts_close(in);
  hts_close(out);

  return 0;
}
