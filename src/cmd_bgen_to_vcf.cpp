#include "cramore.h"
#include "bcf_ordered_writer.h"
#include "tsv_reader.h"
#include "genome_interval.h"
#include <zlib.h>
#include <errno.h>

int32_t cmdBgenToVcf(int32_t argc, char** argv) {
  std::string bgen;
  std::string samplef;
  int32_t verbose = 1000;
  std::string outf;
  std::string region;
  int32_t ifrom = 0;
  int32_t nsnp = INT_MAX;
  bool chrPrefix = false;

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Commands to handle", NULL)
    LONG_STRING_PARAM("bgen",&bgen, "Input VCF/BCF files")
    LONG_STRING_PARAM("sample",&samplef, "Sample ID files")

    LONG_PARAM_GROUP("Region-based access options", NULL)
    LONG_STRING_PARAM("region", &region, "Subset of region to write [CHR]:[BEG]-[END]")

    LONG_PARAM_GROUP("SNP count based access options", NULL)    
    LONG_INT_PARAM("from", &ifrom, "Number of SNPs to skip (Cannot be used with --region)")
    LONG_INT_PARAM("nsnp", &nsnp,  "Maximum number of SNPs to write (Cannot be used with --region)")    

    LONG_PARAM_GROUP("Output Options", NULL)    
    LONG_STRING_PARAM("out",&outf, "Output VCF file")
    LONG_PARAM("chr-prefix",&chrPrefix, "Use chr prefix when representing chromosome names")        
    LONG_INT_PARAM("verbose",&verbose,"Frequency of verbose output (1/n)")    
  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();
  
  // sanity check of input arguments
  if ( bgen.empty() ) {
    error("[E:%s:%d %s] --in, --out are required parameters",__FILE__,__LINE__,__FUNCTION__);
  }

  std::vector<std::string> sample_ids;

  if ( !samplef.empty() ) {
    tsv_reader tr(samplef.c_str());
    int32_t nfields = 0;
    tr.read_line();
    tr.read_line(); // skip first two lines
    while( (nfields = tr.read_line()) > 0 ) {
      sample_ids.push_back(tr.str_field_at(1));
    }
  }
  else {
    //error("--sample is required for v1.1 format");
  }

  GenomeInterval interval;
  if ( !region.empty() ) {
    interval.set(region);
    if ( ( ifrom > 0 ) || ( nsnp < INT_MAX ) )
      error("--from and --nsnp argument cannot be used with --region option togetehr");
  }

  htsFile* wf = hts_open(outf.c_str(), "wz");
  hprintf(wf, "##fileformat=VCFv4.2\n");
  for(int32_t i=1; i <= 22; ++i) {
    hprintf(wf, "##contig=<ID=%s%d>\n", chrPrefix ? "chr" : "", i);
  }
  hprintf(wf, "##contig=<ID=%sX>\n",chrPrefix ? "chr" : "");  
  hprintf(wf, "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Alternate Allele Count\">\n");
  hprintf(wf, "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total Allele Count\">\n");
  hprintf(wf, "##INFO=<ID=GC,Number=G,Type=Float,Description=\"Genotype Count\">\n");    
  hprintf(wf, "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n");  
  hprintf(wf, "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of individuals with non-missing genotypes\">\n");      
  hprintf(wf, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
  hprintf(wf, "##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"Genotype Dosages\">\n");

  FILE* fp = fopen(bgen.c_str(), "r");
  //int32_t int4;
  uint32_t uint4, ret;
  //int32_t  int4;
  uint8_t  uint1;
  //int16_t int2;
  //uint16_t uint2;

  if ( fp == NULL ) error("Cannot open file %s", bgen.c_str());

  if ( ( ret = fread(&uint4, sizeof(uint32_t), 1, fp) ) < 0 ) error("[E:%s:%d %s]",__FILE__,__LINE__,__FUNCTION__);
  int32_t offset = (int32_t) uint4;
  int32_t toskip = offset;
  notice("Offset : %u", uint4);

  if ( ( ret = fread(&uint4, sizeof(uint32_t), 1, fp) ) < 0 ) error("[E:%s:%d %s]",__FILE__,__LINE__,__FUNCTION__);
  int32_t hdrLength = (int32_t)uint4;
  toskip -= sizeof(uint32_t);
  notice("Header block length : %u", hdrLength );
  
  if ( ( ret = fread(&uint4, sizeof(uint32_t), 1, fp) ) < 0 ) error("[E:%s:%d %s]",__FILE__,__LINE__,__FUNCTION__);
  int32_t nvar = (int32_t)uint4;
  toskip -= sizeof(uint32_t);  
  notice("# of Variants : %u", nvar);

  int32_t nsample;

  if ( ( ret = fread(&nsample, sizeof(uint32_t), 1, fp) ) < 0 ) error("[E:%s:%d %s]",__FILE__,__LINE__,__FUNCTION__);
  toskip -= sizeof(uint32_t);  
  notice("# of Samples : %u", nsample);
  
  // 4 byte MAGIC number bgen
  char magicNumber[5] = {0,0,0,0,0};
  if ( ( ret = fread(magicNumber, sizeof(char), 4, fp) ) < 0 ) error("[E:%s:%d %s]",__FILE__,__LINE__,__FUNCTION__);
  toskip -= sizeof(uint32_t);  
  notice("Magic number : %s", magicNumber);
  
  // L_H-20 free data
  if ( hdrLength > 20 ) {
    char* freeData = (char*)malloc(sizeof(char)*(hdrLength-20));
    if ( ( ret = fread(freeData, sizeof(char), hdrLength-20, fp) ) < 0 ) error("[E:%s:%d %s]",__FILE__,__LINE__,__FUNCTION__);
    toskip -= (hdrLength-20);
    notice("Free Data : %d, %s", strlen(freeData), freeData);
  }
  
  // 4 byte set of flags
  if ( ( ret = fread(&uint4, sizeof(uint32_t), 1, fp) ) < 0 ) error("[E:%s:%d %s]",__FILE__,__LINE__,__FUNCTION__);
  toskip -= sizeof(uint32_t);    
  //notice("Compressed SNP blocks : %d", (int32_t)((uint4 >> 30) & 0x03));
  //notice("Layout : %d", (int32_t)((uint4 >> 26) & 0x0f));
  //notice("Sample identifiers : %d", (int32_t)(uint4 & 0x01));
  notice("flags : %x", uint4);
  notice("Compressed SNP blocks : %d", (int32_t)(uint4 & 0x03));
  notice("Layout : %d", (int32_t)((uint4 >> 2) & 0x0f));
  notice("Sample identifiers : %d", (int32_t)((uint4 >> 31) & 0x01));

  if ( ( uint4 & 0x03 ) != 1 )             error("Gzip Compression is not used");
  uint8_t layout = ( uint4 >> 2 ) & 0x0f;
  if ( ( layout != 1 ) && ( layout !=2 ) ) error("Unknown layout is used");

  if ( ( uint4 >> 31 ) & 0x01 )
    error("The routine to read sample identifiers is not implemented yet. You need to specify sample IDs as --sample");
  else if ( sample_ids.empty() )
    warning("Sample ID is not provided using --sample parameter. Assinging IDs as ID_000000000, ....");
  else if ( nsample != (int32_t) sample_ids.size() )
    error("The number of sample IDs provided by --sample parameter (%u) does not match with BGEN file (%d)", sample_ids.size(), nsample);

  while( toskip > 0 ) {
    if ( (ret = fread(&uint1, sizeof(char), 1, fp)) < 0 )
      error("[E:%s:%d %s]",__FILE__,__LINE__,__FUNCTION__);
    --toskip;
  }

  hprintf(wf, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
  
  for(int32_t i=0; i < nsample; ++i) {
    if ( sample_ids.empty() ) {
      hprintf(wf,"\tID_%08d", i+1);      
    }
    else {
      hprintf(wf,"\t%s", sample_ids[i].c_str());
    }
  }
  hprintf(wf, "\n");  
    
  //if ( (int32 = fseek(fp, (int32_t)(offset + 4), SEEK_SET)) < 0 )
  //error("[E:%s:%d %s] Errno is %d. Offset = %u",__FILE__,__LINE__,__FUNCTION__, errno, offset);

  //if ( ret < 0 )
  //error("fseek %d, %d", (int32_t)(offset+4), ret);
  
  Bytef*    comp = (Bytef*) malloc(6 * nsample + 6); //10000000); //char[6 * sample_ids.size()];
  uint32_t  compsz = 6 * nsample + 6;
  uint16_t* gps = new uint16_t[ 3 * nsample ];
  //uint32_t  gpssz = 6 * nsample;
  float*    dss = new float[ nsample ];
  uint8_t*  gts = new uint8_t[ nsample ];
  Bytef*    uncomp = NULL;
  uint32_t  uncompsz = 0;

  int32_t nwritten = 0;
  for(int32_t i=0; i < nvar; ++i) {
    if ( layout == 1 ) {
      int32_t ns;
      if ( ( ret = fread(&ns, sizeof(uint32_t), 1, fp) ) < 0 ) error("[E:%s:%d %s]",__FILE__,__LINE__,__FUNCTION__);
      
      if ( ns != nsample ) {
	error("nsample %d != ns %d", nsample, ns );
      }
    }
    
    int16_t l_varid, l_rsid, l_chr;
    if ( ( ret = fread(&l_varid, sizeof(uint16_t), 1, fp) ) < 0 ) error("[E:%s:%d %s]",__FILE__,__LINE__,__FUNCTION__);
    
    char* varid = new char[l_varid+1];
    if ( l_varid > 0 ) 
      if ( ( ret = fread(varid, sizeof(char), l_varid, fp) ) < 0 ) error("[E:%s:%d %s]",__FILE__,__LINE__,__FUNCTION__);
    varid[l_varid] = '\0';

    //notice("l_varid = %d, varid = %s", l_varid, varid);

    if ( ( ret = fread(&l_rsid, sizeof(uint16_t), 1, fp) ) < 0 ) error("[E:%s:%d %s]",__FILE__,__LINE__,__FUNCTION__);
    char* rsid = new char[l_rsid+1];
    if ( l_rsid > 0 )    
      if ( ( ret = fread(rsid, sizeof(char), l_rsid, fp) ) < 0 ) error("[E:%s:%d %s]",__FILE__,__LINE__,__FUNCTION__);
    rsid[l_rsid] = '\0';

    if ( ( ret = fread(&l_chr, sizeof(uint16_t), 1, fp) ) < 0 ) error("[E:%s:%d %s]",__FILE__,__LINE__,__FUNCTION__);
    char* chr = new char[l_chr+1];
    if ( l_chr > 0 )
      if ( ( ret = fread(chr, sizeof(char), l_chr, fp) ) < 0 ) error("[E:%s:%d %s]",__FILE__,__LINE__,__FUNCTION__);
    chr[l_chr] = '\0';

    // if chr starts with 0..
    if ( chr[0] == '0' ) {
      int32_t nchr = atoi(chr);
      sprintf(chr, "%d", nchr);
    }

    uint32_t pos1;
    if ( ( ret = fread(&pos1, sizeof(uint32_t), 1, fp) ) < 0 ) error("[E:%s:%d %s]",__FILE__,__LINE__,__FUNCTION__);
    uint16_t nalleles = 2;

    if ( layout == 2 ) {
      if ( ( ret = fread(&nalleles, sizeof(uint16_t), 1, fp) ) < 0 ) error("[E:%s:%d %s]",__FILE__,__LINE__,__FUNCTION__); // for 1.1 format
    }

    char buf[65536];
    std::vector<std::string> alleles;
    for(int32_t j=0; j < nalleles; ++j) {
      uint32_t l_allele;
      if ( ( ret = fread(&l_allele, sizeof(uint32_t), 1, fp) ) < 0 ) error("[E:%s:%d %s]",__FILE__,__LINE__,__FUNCTION__);
      //notice("j = %d, l_allele = %u", j, l_allele);
      if ( ( ret = fread(buf, sizeof(char), l_allele, fp) ) < 0 ) error("[E:%s:%d %s]",__FILE__,__LINE__,__FUNCTION__);
      std::string s(buf, l_allele);
      alleles.push_back(s);
    }


    //notice("Lengths : %u %u %u", l_varid, l_rsid, l_chr);
    //notice("Variant ID     : %s", varid);
    //notice("rsID           : %s", rsid);
    //notice("Chromosome     : %s", chr);
    //notice("# Alleles      : %d", nalleles);
    //for(int32_t j=0; j < nalleles; ++j) {
    //  notice(" Allele %d     : %s", j+1, alleles[j].c_str());
    //}

    
    // for layout 1
    // if compressed

    //notice("a/b allele is %s/%s %d", alleles[0].c_str(), alleles[1].c_str(), sizeof(Bytef));

    uint32_t sz;
    if ( ( ret = fread(&sz, sizeof(uint32_t), 1, fp) ) < 0 ) error("[E:%s:%d %s]",__FILE__,__LINE__,__FUNCTION__);
    //notice("Compressed genotype probability size is %u", sz);

    if ( compsz < sz ) {
      compsz = sz + 6;
      comp = (Bytef*)realloc(comp, compsz);
    }

    uint32_t dsz = 6*nsample;
    if ( layout == 2 ) {
      if ( ( ret = fread(&dsz, sizeof(uint32_t), 1, fp) ) < 0 ) error("[E:%s:%d %s]",__FILE__,__LINE__,__FUNCTION__);      
    }
    
    //notice("compressed size   : %u",sz);
    //notice("Uncompressed size : %u",dsz);

    // read compressed data
    if ( ( uint4 = fread(comp, sizeof(Bytef), sz - (layout == 2 ? 4 : 0), fp) ) < 0 ) error("[E:%s:%d %s]",__FILE__,__LINE__,__FUNCTION__);

    if ( i % verbose == 0 ) 
      notice("Processing %d/%d variants, writing %d variants at %s:%u", i+1, nvar, nwritten, chr, pos1);    

    // skip if don't need to uncompress
    if ( !region.empty() ) {
      if ( interval.seq.compare(chr) != 0 ) continue;
      if ( pos1 < (uint32_t)interval.start1 ) continue;
      if ( pos1 > (uint32_t)interval.end1 ) break;
    }
    else if ( i < ifrom ) continue;
    else if ( ( nsnp < INT_MAX) && ( i >= ifrom + nsnp ) ) {
      break;
    }

    // or write the variant

    ++nwritten;

    uLongf usz = dsz;
    uLongf csz = sz;
    
    //Bytef* dst = (Bytef*)malloc(10000000000);

    //notice("Compressed values starts with %x %x %x and end with %x %x %x %x %x", comp[0], comp[1], comp[2], comp[csz-5], comp[csz-4], comp[csz-3], comp[csz-2], comp[csz-1]);

    if ( layout == 2 ) {
      if ( uncompsz < dsz ) {
        uncompsz = dsz;
	uncomp = (Bytef*)realloc(uncomp, dsz);
      }
      int32_t uncompressValue = uncompress(uncomp, &usz, comp, csz);
      if ( uncompressValue != Z_OK )
	error("Error in uncompressing the data %d, %u, %u", uncompressValue, sz, uint4);

      if ( (uint32_t)nsample != *(uint32_t*)(uncomp) )    error("Number of individual in genotype block for %s (%s, %d-th) does not match to %d", varid, rsid, i+1, nsample);
      if ( *(uint16_t*)(uncomp+4) != 2 )        error("Cannot handle multi-allelic variants");
      if ( *(uint8_t*)(uncomp+6) != 2 )         error("Non-diploid cannot be handled");
      if ( *(uint8_t*)(uncomp+7) != 2 )         error("Non-diploid cannot be handled");
      if ( *(uint8_t*)(uncomp+8+nsample) != 0 ) error("Phased data cannot be handled");
      if ( *(uint8_t*)(uncomp+9+nsample) != 8 ) error("non 8-bit probability cannot be handled");

      int32_t an = 0;
      int32_t ac = 0;
      for(int32_t j=0; j < nsample; ++j) {
	if ( uncomp[8+j] == 2 ) {
	  uint8_t* gp = uncomp + 10 + nsample + j*2;
	  uint8_t gp2 = 255-gp[0]-gp[1];
	  an += 2;
	  dss[j] = (gp[1] + gp2 + gp2)/255.0;
	  gts[j] = (gp[0] > gp2) ? (gp[0] > gp[1] ? 1 : 2) : (gp2 > gp[1] ? 3 : 2);	  
	  ac += (gts[j] - 1); 
	}
	else gts[j] = 0;
      }

      hprintf(wf, "%s\t%u\t%s\t%s\t%s\t100\tPASS\tAC=%d;AN=%d;AF=%.3lg;NS=%d\tGT:DS", chr, pos1, rsid, alleles[0].c_str(), alleles[1].c_str(), ac, an, (double)ac/(double)an, (int32_t)(an/2));
      for(int32_t j=0; j < nsample; ++j) {
	if ( gts[j] == 0 ) {
	  hprintf(wf, "\t./.:%.3lg", ac/an*2.0);
	}
	else {
	  hprintf(wf, "\t%d/%d:%.3lg", gts[j] == 3 ? 1 : 0, gts[j] > 1 ? 1 : 0, dss[j]);
	}
      }
      hprintf(wf, "\n");

      //hts_close(wf);

      //if ( i > 100 ) error("foo");
    }
    else if ( layout == 1 ) {
      int32_t uncompressValue = uncompress((Bytef*)gps, &usz, comp, csz);
      //notice("usz = %u", usz);
      
      if ( uncompressValue != Z_OK )
	error("Error in uncompressing the data %d, %u, %u", uncompressValue, sz, uint4);

      
      //notice("Uncompressed genotype probability size is %u", usz);

      //if ( fread(gps, sizeof(uint16_t), nsample*3, fp) < 0 ) error("[E:%s:%d %s]",__FILE__,__LINE__,__FUNCTION__);
      
      // calculate frequencies
      //memset(gps, 0, 3*nsample*sizeof(uint16_t));
      
      double ac = 0;
      int32_t an = 0;
      std::vector<int32_t> imiss;
      for(int32_t j=0; j < nsample; ++j) {
	if ( ( gps[3*j] == 0 ) && ( gps[3*j+1] == 0 ) && ( gps[3*j+2] == 0 ) ) {
	  imiss.push_back(j);
	  gts[j] = 0;
	}
	else {
	  dss[j] = (gps[3*j+1] + 2.0*gps[3*j+2]) / (32767.0);
	  ac += dss[j];
	  an += 2;
	  if ( gps[3*j] > gps[3*j+1] ) {
	    gts[j] = ( gps[3*j] > gps[3*j+2] ? 1 : 3 );
	  }
	  else if ( gps[3*j+1] > gps[3*j+2] ) {
	    gts[j] = 2;
	  }
	  else {
	    gts[j] = 3;
	  }
	}
      }

      //notice("ac=%lf, an=%d",ac,an);
      hprintf(wf, "%s\t%u\t%s\t%s\t%s\t100\tPASS\tAC=%d;AN=%d;AF=%.6lf;NS=%d\tGT:DS", chr, pos1, rsid, alleles[0].c_str(), alleles[1].c_str(), (int32_t)floor(ac+0.5), an, ac/an, an/2);
      for(int32_t j=0; j < nsample; ++j) {
	if ( gts[j] == 0 ) {
	  hprintf(wf, "\t./.:%.3lg", ac/an*2.0);
	}
	else {
	  hprintf(wf, "\t%d/%d:%.3lg", gts[j] == 3 ? 1 : 0, gts[j] > 1 ? 1 : 0, dss[j]);
	}
      }
      hprintf(wf, "\n");
    }

    delete [] varid;
    delete [] rsid;
    delete [] chr;
  }

  hts_close(wf);  

  notice("Finished processing all %d variants, writing %d variants", nvar, nwritten);
  
  delete[] dss;
  delete[] gps;
  delete[] gts;
  if ( comp != NULL ) free(comp);
  if ( uncomp != NULL ) free(uncomp);  

  notice("Analysis finished");

  return 0;
}

