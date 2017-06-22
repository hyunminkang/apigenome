#include "cramore.h"
#include "nuclear_pedigree.h"
#include "bcf_ordered_reader.h"

int32_t getPersonGenoDepth( int32_t* gts, int32_t* dps, NuclearFamilyPerson* pPerson, std::vector<int>& genos, std::vector<int>& depths) {
  genos.clear();
  depths.clear();
  if ( pPerson == NULL ) return 0;
  else {
    int32_t nSamples = 0;
    for(int32_t i=0; i < (int32_t)pPerson->samples.size(); ++i) {
      int32_t idx = pPerson->samples[i]->index;
      if ( idx >= 0 ) {
        int32_t g1 = gts[2*idx];
        int32_t g2 = gts[2*idx+1];
        int32_t geno;
        if ( bcf_gt_is_missing(g1) || bcf_gt_is_missing(g2) ) {
          geno = 0;
        }
        else {
          geno = bcf_alleles2gt(bcf_gt_allele(g1),bcf_gt_allele(g2))+1;
        }
        //fprintf(stderr,"%s %d %d %d %d\n",pPerson->samples[i]->sampleID.c_str(), idx, g1, g2, geno);
        genos.push_back(geno);
        depths.push_back(dps == NULL ? 0 : dps[idx]);
        ++nSamples;
      }
    }
    return nSamples;
  }
}

int32_t cmdVcfMendelDupConc(int32_t argc, char** argv) {
  std::string inVcf;
  std::string inPed;
  std::string region;
  std::string outf;
  int32_t minDP = 0;
  int32_t verbose = 1000;

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input Files", NULL)
    LONG_STRING_PARAM("vcf",&inVcf, "Input VCF/BCF file")
    LONG_STRING_PARAM("ped",&inPed, "Input PED file")    
    LONG_STRING_PARAM("region",&region,"Genomic region to focus on")
    LONG_INT_PARAM("minDP",&minDP,"Minimum genotype depth")
    LONG_INT_PARAM("verbose",&verbose,"Frequency of verbose output (randomly 1/n)")    

    LONG_PARAM_GROUP("Output Files", NULL)
    LONG_STRING_PARAM("out",&outf, "Output file prefix")
  END_LONG_PARAMS();
  
  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  notice("Analysis Started");    
  
  // sanity check of input arguments
  if ( inPed.empty() || outf.empty() || inVcf.empty() ) {
    error("[E:%s:%d %s] --vcf, --out, --ped are required parameters",__FILE__,__LINE__,__FUNCTION__);
  }

  notice("Loading pedigree file %s",inPed.c_str());
  NuclearPedigree* ped = new NuclearPedigree(inPed.c_str());
    
  std::vector<GenomeInterval> intervals;
  if ( !region.empty() ) {
    parse_intervals(intervals, "", region);
  }
  BCFOrderedReader odr(inVcf, intervals);

  int32_t nsamples = bcf_hdr_nsamples(odr.hdr);
  for(int32_t i=0; i < nsamples; ++i)
    ped->setSampleIndex(odr.hdr->samples[i], i);

  notice("In the Original Pedigree: %d families, %d samples, %d unique individuals", (int32_t)ped->famIDmap.size(), (int32_t)ped->smIDmap.size(), ped->numPeople());

  int32_t nremoved = ped->removeSamplesWithoutIndex();
  notice("Removed %d samples not in the VCF file from the pedigree",nremoved);

  notice("Overlapping Samples Only : %d families, %d individuals, %d sequenced samples", (int32_t)ped->famIDmap.size(), ped->numPeople(), ped->numSamplesWithIndex());

  bcf1_t* iv = bcf_init();
  int32_t nread = 0, nskip = 0;

  //int32_t ns = ped->numSamplesWithIndex();
  int32_t* p_gt = NULL;
  int32_t* p_dp = NULL;
  int32_t* p_od = NULL;  
  int32_t np_gt, np_dp, np_od;

  std::vector<int32_t> dadGTs;
  std::vector<int32_t> momGTs;
  std::vector<int32_t> nKids;
  std::vector< std::vector<int32_t> > kidGTs;
  std::vector<int32_t> dadDPs;
  std::vector<int32_t> momDPs;
  std::vector< std::vector<int32_t> > kidDPs;
  std::vector<int32_t> famGTs;
  std::vector<int32_t> trioGTs(3);

  // iterate over each family
  // for each family member, we collect the following metrics
  // (4 x 4) genotype concordance matrix for trios
  // 4^{# dups} matrix for dups
  std::map<std::string, NuclearFamily*>::iterator itF;  
  std::map<NuclearFamily*, FamilyConcordance> famConc;
  std::map<NuclearFamilyPerson*, DupConcordance>    dupConc;
  // print out variant level summary
  std::vector<int32_t> c64, c16;  
  
  int32_t i, j, k;

  htsFile* wf_vf = hts_open((outf+".var.fam.conc").c_str(),"w");
  htsFile* wf_vd = hts_open((outf+".var.dup.conc").c_str(),"w");
  htsFile* wf_if = hts_open((outf+".ind.fam.conc").c_str(),"w");
  htsFile* wf_id = hts_open((outf+".ind.dup.conc").c_str(),"w");

  hprintf(wf_vf,"CHROM\tPOS\tREF\tALT\tTOTAL");
  hprintf(wf_if,"DAD\tMOM\tKID\tTOTAL");
  for(i=0; i < 64; ++i) {
    hprintf(wf_vf,"\tN%d%d%d",(int32_t)(i/16), (int32_t)((i/4) % 4), (int32_t)(i % 4));
    hprintf(wf_if,"\tN%d%d%d",(int32_t)(i/16), (int32_t)((i/4) % 4), (int32_t)(i % 4));    
  }
  hprintf(wf_vd,"CHROM\tPOS\tREF\tALT\tTOTAL");
  hprintf(wf_id,"ID1\tID2\tTOTAL");
  for(i=0; i < 16; ++i) {
    hprintf(wf_vd,"\tN%d%d",(int32_t)(i/4), (int32_t)(i%4));
    hprintf(wf_id,"\tN%d%d",(int32_t)(i/4), (int32_t)(i%4));
  }
  hprintf(wf_vf,"\n");
  hprintf(wf_if,"\n");
  hprintf(wf_vd,"\n");
  hprintf(wf_id,"\n");  
  
  for(itF = ped->famIDmap.begin(); itF != ped->famIDmap.end(); ++itF, ++i) {
    NuclearFamily* pFam = itF->second;

    if ( ( pFam->pDad ) && ( pFam->pDad->samples.size() > 1 ) ) 
      dupConc.insert( std::make_pair(pFam->pDad, DupConcordance((int32_t)pFam->pDad->samples.size())) );
    if ( ( pFam->pMom ) && ( pFam->pMom->samples.size() > 1 ) ) 
      dupConc.insert( std::make_pair(pFam->pMom, DupConcordance((int32_t)pFam->pMom->samples.size())) );
    
    if ( !pFam->pKids.empty() ) {
      famConc.insert( std::make_pair(pFam, FamilyConcordance((int32_t)pFam->pKids.size())) );

      for(i=0; i < (int32_t)pFam->pKids.size(); ++i) {
	if ( pFam->pKids[i]->samples.size() > 1 )
	  dupConc.insert( std::make_pair(pFam->pKids[i], DupConcordance((int32_t)pFam->pKids[i]->samples.size())) );
      }
    }
  }
  
  for(nread=0; odr.read(iv); ++nread) {
    bool skip = false;
    bcf_unpack(iv, BCF_UN_ALL);

    if ( iv->n_allele > 2 )
      skip = true;
    else if ( ( !intervals.empty() ) && ( ( iv->pos + 1 < intervals[0].start1 ) || ( iv->pos + 1 > intervals[0].end1 ) ) )
      skip = true;
    else {
      bool is_vntr = false;
      for(i=0; i < iv->n_allele; ++i) {
	if ( strcmp(iv->d.allele[i],"<VNTR>") == 0 )
	  is_vntr = true;
      }
      if ( is_vntr ) skip = true;      
    }
    if ( skip ) {
      ++nskip;
      continue;
    }

    if ( iv->pos % verbose == 0 ) {
      notice("Reporting whenever the position is a multiple of %d - currently processing [%s %d %s %s]",verbose, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1, iv->d.allele[0], iv->d.allele[1]);    
    }

    // extract GT and AD/DP field
    if ( bcf_get_genotypes(odr.hdr, iv, &p_gt, &np_gt) < 0 )
      error("[E:%s:%d %s] Cannot parse GT field",__FILE__,__LINE__,__FUNCTION__);

    if ( minDP > 0 ) {
      if ( bcf_get_format_int32(odr.hdr, iv, "DP", &p_dp, &np_dp) < 0 ) {
	if ( bcf_get_format_int32(odr.hdr, iv, "AD", &p_dp, &np_dp) < 0 ) {
	  error("[E:%s:%d %s] Cannot parse AD or DP field",__FILE__,__LINE__,__FUNCTION__);	
	}
	else if ( bcf_get_format_int32(odr.hdr, iv, "OD", &p_od, &np_od) < 0 ) {
	  error("[E:%s:%d %s] Cannot parse AD or DP field",__FILE__,__LINE__,__FUNCTION__);		  
	}
	
	// if AD field is available, use their sum as depth (assuming biallelics);
	for(i=0; i < nsamples; ++i) {
	  p_dp[i] = p_dp[2*i] + p_dp[2*i+1] + p_od[i];
	}
      }
    }

    FamilyConcordance vFam(1);
    DupConcordance vDup(2);
    
    for(itF = ped->famIDmap.begin(); itF != ped->famIDmap.end(); ++itF, ++i) {
      NuclearFamily* pFam = itF->second;

      // calculate genotype concordance
      // first get genotypes
      famGTs.clear();
      
      int32_t nDad = getPersonGenoDepth( p_gt, p_dp, pFam->pDad, dadGTs, dadDPs);
      for(i=0; i < nDad; ++i) {
	if ( dadDPs[i] < minDP )
	  dadGTs[i] = 0;
      }
      famGTs.push_back(nDad > 0 ? dadGTs[0] : 0);
      trioGTs[0] = (nDad > 0 ? dadGTs[0] : 0);

      int32_t nMom = getPersonGenoDepth( p_gt, p_dp, pFam->pMom, momGTs, momDPs);
      for(i=0; i < nMom; ++i) {
	if ( momDPs[i] < minDP )
	  momGTs[i] = 0;
      }
      famGTs.push_back(nMom > 0 ? momGTs[0] : 0);
      trioGTs[1] = (nMom > 0 ? momGTs[0] : 0);      

      if ( nDad > 1 ) {
	dupConc[pFam->pDad].addGenotype(dadGTs);
	for( i=1; i < nDad; ++i )
	  for( k=0; k < i; ++k) 
	    vDup.addGenotype(dadGTs[k],dadGTs[i]);
      }
      if ( nMom > 1 ) {
	dupConc[pFam->pMom].addGenotype(momGTs);
	for( i=1; i < nMom; ++i )
	  for( k=0; k < i; ++k) 
	    vDup.addGenotype(momGTs[k],momGTs[i]);	
      }

      if ( !pFam->pKids.empty() ) {
	nKids.resize(pFam->pKids.size());
	kidGTs.resize(pFam->pKids.size());
	kidDPs.resize(pFam->pKids.size());
	for(j=0; j < (int32_t)pFam->pKids.size(); ++j) {
	  nKids[j] = getPersonGenoDepth( p_gt, p_dp, pFam->pKids[j], kidGTs[j], kidDPs[j]);
	  for(i=0; i < nKids[j]; ++i) {
	    if ( kidDPs[j][i] < minDP )
	      kidGTs[j][i] = 0;
	    trioGTs[2] = kidGTs[j][i];
	    vFam.addGenotype(trioGTs);
	  }
	  famGTs.push_back( (nKids[j] > 0) ? kidGTs[j][0] : 0);
	}
	
	for(j=0; j < (int32_t)pFam->pKids.size(); ++j) {
	  if ( nKids[j] > 1 ) {
	    dupConc[pFam->pKids[j]].addGenotype(kidGTs[j]);
	    for( i=1; i < nKids[j]; ++i )
	      for( k=0; k < i; ++k) 
		vDup.addGenotype(kidGTs[j][k],kidGTs[j][i]);
	  }
	}
	
	// get the duplicate concordance and trio concordance
	famConc[pFam].addGenotype(famGTs);
      }
    }

    std::string hdr;
    int32_t total = vFam.fillTrioCount(0,c64);
    catprintf(hdr, "%s\t%d\t%s\t%s",bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1, iv->d.allele[0], iv->d.allele[1]);    
    if ( total > 0 ) {
      printTrioDupCount(wf_vf, hdr, c64);
    }

    total = vDup.fillDupCount(0,1,c16);
    if ( total > 0 ) {    
      printTrioDupCount(wf_vd, hdr, c16);
    }
  }

  for(itF = ped->famIDmap.begin(); itF != ped->famIDmap.end(); ++itF, ++i) {
    NuclearFamily* pFam = itF->second;
    FamilyConcordance& iFam = famConc[pFam];
    int32_t total;

    std::string dadID = pFam->pDad ? pFam->pDad->samples[0]->sampleID : ".";
    std::string momID = pFam->pMom ? pFam->pMom->samples[0]->sampleID : ".";

    std::string hdr;     
    // print duplicates for dad
    if ( ( pFam->pDad ) && ( pFam->pDad->samples.size() > 1 ) ) {
      DupConcordance& iDup = dupConc[pFam->pDad];
      for(i=1; i < (int32_t)pFam->pDad->samples.size(); ++i) {
	for(k=0; k < i; ++k) {
	  total = iDup.fillDupCount(k,i,c16);
	  if ( total > 0 ) {
	    hdr.clear();
	    catprintf(hdr,"%s\t%s",pFam->pDad->samples[k]->sampleID.c_str(),pFam->pDad->samples[i]->sampleID.c_str());
	    printTrioDupCount(wf_id, hdr, c16);
	  }
	}
      }
    }

    // print duplicates for mom
    if ( ( pFam->pMom ) && ( pFam->pMom->samples.size() > 1 ) ) {
      DupConcordance& iDup = dupConc[pFam->pMom];
      for(i=1; i < (int32_t)pFam->pMom->samples.size(); ++i) {
	for(k=0; k < i; ++k) {
	  total = iDup.fillDupCount(k,i,c16);
	  if ( total > 0 ) {
	    hdr.clear();	    
	    catprintf(hdr,"%s\t%s",pFam->pMom->samples[k]->sampleID.c_str(),pFam->pMom->samples[i]->sampleID.c_str());
	    printTrioDupCount(wf_id, hdr, c16);
	  }
	}
      }
    }    

    if ( !pFam->pKids.empty() ) {
      for(j=0; j < (int32_t)pFam->pKids.size(); ++j) {
	std::string kidID = pFam->pKids[j]->samples[0]->sampleID;
	total = iFam.fillTrioCount(j,c64);
	if ( total > 0 ) {
	  hdr.clear();	  
	  catprintf(hdr,"%s\t%s\t%s",dadID.c_str(), momID.c_str(), kidID.c_str());
	  printTrioDupCount(wf_if, hdr, c64);
	}
	if ( pFam->pKids[j]->samples.size() > 1 ) {
	  DupConcordance& iDup = dupConc[pFam->pKids[j]];
	  for(i=1; i < (int32_t)pFam->pKids[j]->samples.size(); ++i) {
	    for(k=0; k < i; ++k) {
	      total = iDup.fillDupCount(k,i,c16);
	      if ( total > 0 ) {
		hdr.clear();		
		catprintf(hdr,"%s\t%s",pFam->pKids[j]->samples[k]->sampleID.c_str(),pFam->pKids[j]->samples[i]->sampleID.c_str());
		printTrioDupCount(wf_id, hdr, c16);
	      }
	    }
	  }	  
	}
      }
    }
  }

  hts_close(wf_vf);
  hts_close(wf_vd);
  hts_close(wf_if);
  hts_close(wf_id);

  bcf_destroy(iv);

  notice("Analysis Finished");

  return 0;
}

