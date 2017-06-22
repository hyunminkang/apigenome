#include "cramore.h"
#include "bcf_ordered_writer.h"
#include "tsv_reader.h"
#include <zlib.h>

int32_t cmdBgenToVcf(int32_t argc, char** argv) {
  std::string bgen;
  std::string samplef;
  int32_t verbose = 1000;
  std::string outf;

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input VCFs", NULL)
    LONG_STRING_PARAM("bgen",&bgen, "Input VCF/BCF files")
    LONG_STRING_PARAM("sample",&samplef, "Sample ID files")
    LONG_STRING_PARAM("out",&outf, "Output VCF file")    
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
    error("--sample is required for v1.1 format");
  }

  htsFile* wf = hts_open(outf.c_str(), "w");
  hprintf(wf, "##fileformat=VCFv4.2\n");
  for(int32_t i=1; i <= 22; ++i) {
    hprintf(wf, "##contig=<ID=%d>\n", i);
  }
  hprintf(wf, "##contig=<ID=X>\n");  
  hprintf(wf, "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Alternate Allele Count\">\n");
  hprintf(wf, "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total Allele Count\">\n");
  hprintf(wf, "##INFO=<ID=GC,Number=G,Type=Float,Description=\"Genotype Count\">\n");    
  hprintf(wf, "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n");  
  hprintf(wf, "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of individuals with non-missing genotypes\">\n");      
  hprintf(wf, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
  hprintf(wf, "##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"Genotype Dosages\">\n");
  hprintf(wf, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
  for(int32_t i=0; i < (int32_t)sample_ids.size(); ++i) {
    hprintf(wf,"\t%s", sample_ids[i].c_str());
  }
  hprintf(wf, "\n");

  FILE* fp = fopen(bgen.c_str(), "r");
  //int32_t int4;
  uint32_t uint4;
  //int16_t int2;
  //uint16_t uint2;

  if ( fp == NULL ) error("Cannot open file %s", bgen.c_str());

  if ( fread(&uint4, sizeof(uint32_t), 1, fp) < 0 ) error("[E:%s:%d %s]",__FILE__,__LINE__,__FUNCTION__);

  int32_t offset = (int32_t) uint4;
  notice("Offset : %u", uint4);

  if ( fread(&uint4, sizeof(uint32_t), 1, fp) < 0 ) error("[E:%s:%d %s]",__FILE__,__LINE__,__FUNCTION__);
  if ( fread(&uint4, sizeof(uint32_t), 1, fp) < 0 ) error("[E:%s:%d %s]",__FILE__,__LINE__,__FUNCTION__);

  int32_t nvar = (int32_t)uint4;
  Bytef*    comp = (Bytef*) malloc(10000000); //char[6 * sample_ids.size()];
  uint16_t* gps = new uint16_t[ 3 * sample_ids.size() ];
  float*    dss = new float[ sample_ids.size() ];
  uint8_t*  gts = new uint8_t[ sample_ids.size() ];     

  if ( fseek(fp, offset + 4, SEEK_SET) < 0 ) error("[E:%s:%d %s]",__FILE__,__LINE__,__FUNCTION__);

  for(int32_t i=0; i < nvar; ++i) {
    int32_t nsample;
    if ( fread(&nsample, sizeof(uint32_t), 1, fp) < 0 ) error("[E:%s:%d %s]",__FILE__,__LINE__,__FUNCTION__);
    if ( nsample != (int32_t)sample_ids.size() ) {
      error("nsample %d != sample_ids.size() = %u", nsample, sample_ids.size());
    }
    
    int32_t l_varid, l_rsid, l_chr;
    if ( fread(&l_varid, sizeof(uint16_t), 1, fp) < 0 ) error("[E:%s:%d %s]",__FILE__,__LINE__,__FUNCTION__);
    
    char* varid = new char[l_varid+1];
    if ( fread(varid, sizeof(char), l_varid, fp) < 0 ) error("[E:%s:%d %s]",__FILE__,__LINE__,__FUNCTION__);
    varid[l_varid] = '\0';

    if ( fread(&l_rsid, sizeof(uint16_t), 1, fp) < 0 ) error("[E:%s:%d %s]",__FILE__,__LINE__,__FUNCTION__);
    char* rsid = new char[l_rsid+1];
    if ( fread(rsid, sizeof(char), l_rsid, fp) < 0 ) error("[E:%s:%d %s]",__FILE__,__LINE__,__FUNCTION__);
    rsid[l_rsid] = '\0';

    if ( fread(&l_chr, sizeof(uint16_t), 1, fp) < 0 ) error("[E:%s:%d %s]",__FILE__,__LINE__,__FUNCTION__);
    char* chr = new char[l_chr+1];
    if ( fread(chr, sizeof(char), l_chr, fp) < 0 ) error("[E:%s:%d %s]",__FILE__,__LINE__,__FUNCTION__);
    chr[l_chr] = '\0';

    uint32_t pos1;
    if ( fread(&pos1, sizeof(uint32_t), 1, fp) < 0 ) error("[E:%s:%d %s]",__FILE__,__LINE__,__FUNCTION__);
    uint16_t nalleles = 2;
    //if ( fread(&nalleles, sizeof(uint16_t), 1, fp) < 0 ) error("[E:%s:%d %s]",__FILE__,__LINE__,__FUNCTION__); // for 1.1 format

    notice("%s %u %s %u %d", chr, pos1, rsid, (uint32_t)nalleles, nsample);
    
    char buf[65536];
    std::vector<std::string> alleles;
    for(int32_t j=0; j < nalleles; ++j) {
      uint32_t l_allele;
      if ( fread(&l_allele, sizeof(uint32_t), 1, fp) < 0 ) error("[E:%s:%d %s]",__FILE__,__LINE__,__FUNCTION__);
      notice("j = %d, l_allele = %u", j, l_allele);
      if ( fread(buf, sizeof(char), l_allele, fp) < 0 ) error("[E:%s:%d %s]",__FILE__,__LINE__,__FUNCTION__);
      std::string s(buf, l_allele);
      alleles.push_back(s);
    }

    // for layout 1
    // if compressed

    notice("a/b allele is %s/%s %d", alleles[0].c_str(), alleles[1].c_str(), sizeof(Bytef));

    uint32_t sz;
    if ( fread(&sz, sizeof(uint32_t), 1, fp) < 0 ) error("[E:%s:%d %s]",__FILE__,__LINE__,__FUNCTION__);
    notice("Compressed genotype probability size is %u", sz);

    if ( ( uint4 = fread(comp, sizeof(Bytef), sz, fp) ) < 0 ) error("[E:%s:%d %s]",__FILE__,__LINE__,__FUNCTION__);
    uLongf usz = 0;
    uLongf csz = sz;
    Bytef* dst = (Bytef*)malloc(10000000000);
    int32_t uncompressValue = uncompress(dst, &usz, comp, csz);

    notice("usz = %u", usz);

    if ( uncompressValue != Z_OK )
      error("Error in uncompressing the data %d, %u", uncompressValue, uint4);

    notice("Uncompressed genotype probability size is %u", usz);

    //if ( fread(gps, sizeof(uint16_t), nsample*3, fp) < 0 ) error("[E:%s:%d %s]",__FILE__,__LINE__,__FUNCTION__);

    // calculate frequencies
    double ac = 0;
    int32_t an = 0;
    std::vector<int32_t> imiss;
    for(int32_t j=0; j < nsample; ++j) {
      if ( ( gps[3*j] == 0 ) && ( gps[3*j+1] == 0 ) && ( gps[3*j+2] == 0 ) ) {
	imiss.push_back(j);
	gts[j] = 0;
      }
      else {
	dss[j] = (gps[3*j+1] + 2.0*gps[3*j+2]) / (32768.0);
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

    notice("ac=%lf, an=%d",ac,an);

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

    delete [] varid;
    delete [] rsid;
    delete [] chr;
  }
  
  delete[] dss;
  delete[] gps;
  delete[] gts;

  hts_close(wf);

  notice("Analysis finished");

  return 0;
}

