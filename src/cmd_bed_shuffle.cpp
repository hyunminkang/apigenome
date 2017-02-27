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
#include "reference_sequence.h"

int32_t cmdBedShuffle(int32_t argc, char** argv) {
  std::string inBed;
  std::string outBed;
  std::string maskBed;
  std::string refFasta;
  int32_t maxTries = 1000;
  int32_t verbose = 100000;
  int32_t seed = 0;
  //int32_t gcWindow = 100;
  //int32_t gcBins = 100;
  //std::string gcFile;

  paramList pl;
  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Required parameters", NULL)
    LONG_STRING_PARAM("bed",&inBed,"Input BED file representing the positive labels")
    LONG_STRING_PARAM("out",&outBed,"Output BED file name")
    LONG_STRING_PARAM("ref",&refFasta,"FASTA format reference genome")

    LONG_PARAM_GROUP("Mask parameters", NULL)    
    LONG_STRING_PARAM("mask",&maskBed,"BED file used for masking genomic regions")
    LONG_INT_PARAM("max-tries",&maxTries,"Maximum number of repeated shuffling attempts per locus")
    LONG_INT_PARAM("verbose",&verbose,"Frequency of verbose output")
    LONG_INT_PARAM("seed",&seed,"Seed for random number generator")      
  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  // sanity check of input arguments
  if ( inBed.empty() || outBed.empty() || refFasta.empty() )
    error("[E:%s:%d %s] --pos --neg, --out, --ref are required parameters",__FILE__,__LINE__,__FUNCTION__);

  // read reference sequences
  ReferenceSequence ref(refFasta);
  genomeLoci maskLoci;

  if ( seed == 0 )
    srand(time(NULL));
  else
    srand(seed);

  //std::set<string> maskChrs;

  if ( !maskBed.empty() ) {
    htsFile* hp = hts_open(maskBed.c_str(),"r");
    if ( hp == NULL )
      error("[E:%s:%d %s] Cannot open file %s for reading",__FILE__,__LINE__,__FUNCTION__, maskBed.c_str());
    kstring_t str = {0,0,0};
    int32_t lstr = 0;
    int32_t nfields = 0;
    int32_t* fields = NULL;
    int32_t i;
    // model list is assumed to have [INFO_KEY] [MODEL_FILE] [INFO_DESCRIPTION = INFO_KEY if empty]
    for( i=0; ( lstr = hts_getline(hp, KS_SEP_LINE, &str) ) >= 0; ++i ) {
      fields = ksplit(&str, 0, &nfields);
      if ( nfields < 3 )
	error("[E:%s:%d %s] Less than three columns observed in line %d of %s",__FILE__,__LINE__,__FUNCTION__, i+1, maskBed.c_str());      

      // typically bed files have beg0 and end0
      maskLoci.add(&str.s[fields[0]], atoi(&str.s[fields[1]])+1, atoi(&str.s[fields[2]]));
      //maskChrs.insert(&str.s[fields[0]]);

      if ( i % verbose == 0 )
	notice("Processing %d lines to mask at %s:%d-%d in %s..", i, &str.s[fields[0]], &str.s[fields[1]], &str.s[fields[2]], maskBed.c_str());
    }
    hts_close(hp);
    
    notice("Processed %d lines to mask, maxLength = %d", i, maskLoci.maxLength);
    
    maskLoci.resolveOverlaps();

    notice("After removing overlaps, %d intervals to mask, maxLength = %d", (int32_t)maskLoci.loci.size(), maskLoci.maxLength);    
  }

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
  // model list is assumed to have [INFO_KEY] [MODEL_FILE] [INFO_DESCRIPTION = INFO_KEY if empty]
  for( int32_t i=0; ( lstr = hts_getline(hp, KS_SEP_LINE, &str) ) >= 0; ++i ) {
    fields = ksplit(&str, 0, &nfields);
    if ( nfields < 3 )
      error("[E:%s:%d %s] Less than three columns observed in line %d of %s",__FILE__,__LINE__,__FUNCTION__, i+1, outBed.c_str());
    
    // leave the chromosome, shiffle the
    const char* chr = &str.s[fields[0]];
    int32_t szchr = ref.fetch_seq_len(chr);
    int32_t beg1 = atoi(&str.s[fields[1]])+1;
    int32_t end0 = atoi(&str.s[fields[2]]);
    bool notyet = true;
    int32_t tries = 0;    
    while( notyet ) {
      int32_t newbeg1 = ((uint64_t)beg1 + (uint64_t)rand() * rand()) % ( szchr - (beg1-end0+1) ) + 1;
      int32_t newend0 = end0 + newbeg1 - beg1;
      if ( !maskLoci.overlaps(chr,newbeg1,newend0) ) {}
      //notice("%s:%d-%d does not overlap with mask (original: %s:%d-%d)", chr, newbeg1, newend0, chr, beg1, end0);
      else if ( outLoci.overlaps(chr,newbeg1,newend0) ) {}
        //notice("%s:%d-%d overlaps with output  (original: %s:%d-%d)", chr, newbeg1, newend0, chr, beg1, end0);
      else {
	hprintf(wp, "%s\t%d\t%d\t%s:%d-%d\n", chr, newbeg1, newend0, chr, beg1, end0);
	outLoci.add(chr, newbeg1, newend0);
	sumTries += (tries+1);
	break;
      }
      ++tries;
      if ( tries > maxTries )
	error("More than %d sampling for %s:%d-%d attemped. i=%d", maxTries, chr, beg1, end0, i);
    }

    if ( i % verbose == 0 )
      notice("Processing %d lines to shffule at %s:%d-%d in %s..", i, chr, beg1, end0, inBed.c_str());
  }
  hts_close(hp);
  hts_close(wp);

  notice("Processes %u lines to shffule. Average sampling attempt per locus is %lf", outLoci.loci.size(), sumTries/(double)outLoci.loci.size());

  return 0;
}

