#include "cramore.h"

/* cmd_bed_matched_shuffle.cpp
 *
 * Copyright (C) 2016 Hyun Min Kang
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include "genomeLoci.h"
#include "fastaGC.h"
#include "reference_sequence.h"

// matched shuffling works this way
// 1. Each genomic coordinate is assigned with a locusCovariate
// 2. The center of the region is matched

class locusCovariate {
public:
  int32_t gcBin;
  int32_t fracBin;
  bool    inMask;

  locusCovariate(int32_t _gcBin, int32_t _fracBin, bool _inMask) : gcBin(_gcBin), fracBin(_fracBin), inMask(_inMask) {}

  bool operator==(const locusCovariate& lcov) const {
    return ( ( gcBin == lcov.gcBin ) && ( fracBin == lcov.fracBin ) && ( inMask == lcov.inMask ) );
  }

  bool operator<(const locusCovariate& lcov) const {
    if ( gcBin == lcov.gcBin ) {
      if ( fracBin == lcov.fracBin )
	return inMask < lcov.inMask;
      else
	return fracBin < lcov.fracBin;
    }
    else
      return gcBin < lcov.gcBin;
  }
};

class baseSampler {
public:
  std::vector<genomeLocus> loci;
  std::vector<int32_t>     baseidx;
  
  void add(const genomeLocus& l) {
    int32_t ls = loci.size();
    loci.push_back(l);
    for(int32_t i=l.beg1; i <= l.end0; ++i) {
      baseidx.push_back(ls);
    }
  }

  void sample(std::string& chr, int64_t& pos1) {
    int32_t bi = baseidx[(int32_t)floor((rand()+0.5) / (RAND_MAX+1.0) * (double) baseidx.size())];
    chr = loci[bi].chrom;
    pos1 = (int64_t)loci[bi].beg1 + (int64_t)floor((rand() + 0.5) / (RAND_MAX+1.0) * (loci[bi].end0 - loci[bi].beg1));
  }
};

int32_t cmdBedMatchedShuffle(int32_t argc, char** argv) {
  std::string inBed;
  std::string outBed;
  std::string gcFile;
  std::string repeatBed;
  std::string maskBed;  
  int32_t gcBin = 1;
  double repeatFracBin = 0.01;
  //int32_t maxTries = 1000;
  int32_t verbose = 100000;
  int32_t seed = 0;

  paramList pl;
  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Required parameters", NULL)
    LONG_STRING_PARAM("in",&inBed,"Input BED file representing the positive labels")
    LONG_STRING_PARAM("out",&outBed,"Output BED file name")
    LONG_STRING_PARAM("gc",&gcFile,"GC content find generated from cramore fasta_gc_content")

    LONG_PARAM_GROUP("Options parameters", NULL)    
    LONG_INT_PARAM("bin-gc",&gcBin,"Unit of binning GC content values for matching")
    LONG_STRING_PARAM("repeat",&repeatBed,"Repeat mask BED file")        
    LONG_DOUBLE_PARAM("bin-frac",&repeatFracBin,"Unit of binning repeat mask fraction for matching")
    LONG_STRING_PARAM("mask",&maskBed,"Accessibility mask to focus on")            
    LONG_INT_PARAM("verbose",&verbose,"Frequency of verbose output")
    LONG_INT_PARAM("seed",&seed,"Seed for random number generator")      
  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  // sanity check of input arguments
  if ( inBed.empty() || outBed.empty() || gcFile.empty() )
    error("[E:%s:%d %s] --in --out, --gc are required parameters",__FILE__,__LINE__,__FUNCTION__);

  notice("Loading GC content file from %s", gcFile.c_str());

  // read GC contents profile
  fastaGC fGC;
  fGC.on_memory = true;
  if ( !fGC.openGC(gcFile.c_str()) )
    error("Failed loading GC content file %s", gcFile.c_str());

  notice("Loading BED file %s of repeat mask", repeatBed.c_str());

  genomeLoci maskLoci;
  if ( !maskBed.empty() )
    maskLoci.openBED( maskBed.c_str() );
  
  // read BED file to mask
  genomeLoci repeatLoci;
  if ( !repeatBed.empty() )
    repeatLoci.openBED( repeatBed.c_str() );

  if ( seed == 0 )
    srand(time(NULL));
  else
    srand(seed);

  // try to generate a map between (GCbin, fracBin)
  // memory needed is roughly [# of bins] * (4)
  notice("Generating a map of regions stratified by covariates");
  std::vector<uint32_t> chr_nbins;
  std::vector<uint32_t> cum_nbins;
  for(int32_t i=0; i < (int32_t)fGC.seqnames.size(); ++i) {
    int32_t nbins = (int32_t)ceil((double)fGC.seqlens[i]/(double)fGC.sliding_unit);    
    chr_nbins.push_back( nbins );
    cum_nbins.push_back( i == 0 ? 0 : cum_nbins[i-1] + chr_nbins[i-1] );
  }
  cum_nbins.push_back(cum_nbins.back() + chr_nbins.back());

  //std::map< locusCovariate, std::set<uint32_t> > lcset;
  std::map< locusCovariate, std::vector<uint32_t> > lcvec;  
  std::map<std::string,int32_t> chr2ichr; 
  for(int32_t i=0; i < (int32_t)fGC.seqnames.size(); ++i) {
    const char* chr = fGC.seqnames[i].c_str();
    uint16_t* m   = fGC.mem_gcs[i];
    int32_t nbins = (int32_t)ceil((double)fGC.seqlens[i]/(double)fGC.sliding_unit);
    //int32_t prevpos = -1;
    int32_t gc  = m[0] / gcBin;
    int32_t pos1  = 1 + fGC.sliding_unit / 2;

    // mflag and rflag represents whether the chromosome is in the mask/repeat loci
    bool mflag = maskLoci.hasChrom(chr);
    bool rflag = repeatLoci.hasChrom(chr);

    notice("i=%d, chr=%s, gc = %d, pos1 = %d, nbins = %d, mflag = %d, rflag = %d", i,chr, gc, pos1, nbins, mflag, rflag);        

    if ( ( !maskBed.empty() ) && ( !mflag ) )   continue;
    if ( ( !repeatBed.empty() ) && ( !rflag ) ) continue;

    chr2ichr[chr] = i;

    int32_t frac = 0;
    if ( !repeatLoci.isend() ) {
      int32_t no = repeatLoci.it->overlapBases(chr, pos1-fGC.window_size,pos1+fGC.window_size);
      frac = (int32_t)(no / repeatFracBin / (fGC.window_size*2.0+1.0));
    }

    bool inMask = maskLoci.overlaps(chr, pos1, pos1);
    locusCovariate lcov(gc, frac, inMask);
    lcvec[lcov].push_back(cum_nbins[i]);
    for(int32_t j=1; j < nbins; ++j) {
      pos1 += fGC.sliding_unit;

      //notice("j=%d, pos1=d",j,pos1);      
      if ( j % 10000000 == 0 )
	notice("Processing %s:%d", chr, pos1);      
      
      inMask = maskLoci.overlaps(chr, pos1, pos1);
      gc = m[j]/gcBin;
      frac = repeatLoci.it->overlapBases(chr, pos1-fGC.window_size,pos1+fGC.window_size) / (double)repeatFracBin / (fGC.window_size*2.0+1.0);

      lcov.gcBin = gc;
      lcov.fracBin = frac;
      lcov.inMask = inMask;

      lcvec[lcov].push_back(cum_nbins[i]+j);
    }
  }

  //notice("Finished processing %d loci and %ld bases in total");

  genomeLoci outLoci;
  htsFile* hp = hts_open(inBed.c_str(),"r");
  htsFile* wp = hts_open(outBed.c_str(),"w");
  if ( hp == NULL )
    error("[E:%s:%d %s] Cannot open file %s for reading",__FILE__,__LINE__,__FUNCTION__, inBed.c_str());
  
  if ( wp == NULL )
    error("[E:%s:%d %s] Cannot open file %s for writing",__FILE__,__LINE__,__FUNCTION__, outBed.c_str());  
  
  kstring_t str = {0,0,0};
  int32_t lstr = 0;
  int32_t nfields = 0;
  int32_t* fields = NULL;
  int64_t sumTries = 0;

  for( int32_t i=0; ( lstr = hts_getline(hp, KS_SEP_LINE, &str) ) >= 0; ++i ) {
    fields = ksplit(&str, 0, &nfields);
    if ( nfields < 3 )
      error("[E:%s:%d %s] Less than three columns observed in line %d of %s",__FILE__,__LINE__,__FUNCTION__, i+1, outBed.c_str());
    
    // leave the chromosome, shuffle the
    const char* chr = &str.s[fields[0]];
    //int32_t szchr = ref.fetch_seq_len(chr);
    int32_t beg1 = atoi(&str.s[fields[1]])+1;
    int32_t end0 = atoi(&str.s[fields[2]]);
    //bool notyet = true;
    //int32_t tries = 0;
    std::string s_chr;
    //int64_t s_pos1 = 0;
    int32_t med1 = (beg1+end0)/2;

    if ( chr2ichr.find(chr) == chr2ichr.end() )
      error("Chromosome %s cannot be used due to lack of mask data",chr);
      
    int32_t ichr = chr2ichr[chr];
    uint32_t ubin = cum_nbins[ichr] + med1/fGC.sliding_unit;
    int32_t gc = fGC.mem_gcs[ichr][ubin/gcBin];
    int32_t frac = repeatLoci.it->overlapBases(chr, med1-fGC.window_size,med1+fGC.window_size) / (double)repeatFracBin / (fGC.window_size*2.0+1.0);
    bool inMask = maskLoci.overlaps(chr, med1, med1);
    
    locusCovariate lcov(gc, frac, inMask);

    //std::set<uint32_t>::iterator it = lcset[lcov].begin();
    //std::advance(it, (int)((rand() + 0.5) / (RAND_MAX+1.0) * lcset[lcov].size()));
    int32_t ridx = (int)((rand() + 0.5) / (RAND_MAX+1.0) * lcvec[lcov].size());
    uint32_t rbin = lcvec[lcov][ridx];
    lcvec[lcov][ridx] = lcvec[lcov].back();
    lcvec[lcov].pop_back();

    // convert bin to region again
    //uint32_t rbin = *it;
    for(int32_t j=0; j < (int32_t)fGC.seqnames.size(); ++j) {
      if ( ( cum_nbins[j] <= rbin ) && ( rbin < cum_nbins[j+1] ) ) {
	int32_t med2 = (rbin-cum_nbins[j])*fGC.sliding_unit;
	hprintf(wp, "%s\t%d\t%d\t%s:%d-%d\n", fGC.seqnames[j].c_str(), med2-(end0-beg1)/2, med2-(end0-beg1)/2+(end0-beg1),chr,beg1,end0);
	//lcset[lcov].erase(it);
	break;
      }
    }    
    if ( i % verbose == 0 )
      notice("Processing %d lines to shffule at %s:%d-%d in %s..", i, chr, beg1, end0, inBed.c_str());
  }
  hts_close(hp);
  hts_close(wp);

  notice("Processes %u lines to shuffle. Average sampling attempt per locus is %lf", outLoci.loci.size(), sumTries/(double)outLoci.loci.size());
							 
  return 0;
}

