#include "cramore.h"
#include "sam_filtered_reader.h"
#include "sam_ordered_writer.h"

int32_t cmdCramSimulContam(int32_t argc, char** argv) {
  SAMFilteredReader sr1;
  SAMFilteredReader sr2;
  std::string outPrefix;
  double fracSample1 = 0;
  double fracSample2 = 0;  
  int32_t seed = 0;
  int32_t verbose = 100000;

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Options for input/output SAM/BAM/CRAM", NULL)
    LONG_STRING_PARAM("sam1",&sr1.sam_file_name, "Intended input SAM/BAM/CRAM file to be sampled at (1-alpha) proportion")
    LONG_STRING_PARAM("sam2",&sr2.sam_file_name, "Intended input SAM/BAM/CRAM file to be sampled at (alpha) proportion")
    LONG_STRING_PARAM("out",&outPrefix, "Output SAM/BAM/CRAM file name to write")    

    LONG_PARAM_GROUP("Options for specifying contamination fraction", NULL)
    LONG_DOUBLE_PARAM("frac1",&fracSample1, "Fraction of sampling from the first SAM/BAM/CRAM")
    LONG_DOUBLE_PARAM("frac2",&fracSample2, "Fraction of sampling from the second SAM/BAM/CRAM")    
    LONG_INT_PARAM("seed",&seed,"Randomization seed")

    LONG_PARAM_GROUP("Analysis Options", NULL)
    LONG_INT_PARAM("verbose", &verbose, "Verbosity parameters")
  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  sr1.verbose = verbose;
  sr2.verbose = verbose;   

  sr1.set_buffer_size(1);
  sr1.init_params();
  sr2.set_buffer_size(1);
  sr2.init_params();


  if ( seed == 0 ) seed = (int32_t)time(NULL);

  // scan VCF and CRAM simultaneously
  // read a variant first
  SAMOrderedWriter sow(outPrefix.c_str());
  sow.set_hdr(sr1.hdr);
  //sam_hdr_merge(sow.hdr, sr2.hdr);
  sow.write_hdr();

  /*
  uint16_t humi;
  double dumi;


  int32_t n1 = 0, t1 = 0;
  int32_t n2 = 0, t2 = 0;
  
  while( sr1.read() ) { // read SAM file
    bam1_t* b1 = sr1.cursor();
    bam1_t* b2 = sr2.eof ? NULL : sr2.cursor();
    
    while( ( b2 != NULL ) && ( ( b2->core.tid < b1->core.tid ) || ( ( b2->core.tid == b1->core.tid ) && ( ( b2->pos < b1->pos ) || ( ( b2->pos == b1-pos ) && ( strcmp(bam_name(b2), bam_name(b1)) < 0 ) ) ) ) ) ) {
      humi = (str_hash(bam_name(b2)) % UINT16_MAX);
      dumi = (double)(humi+0.5) / (double)UINT16_MAX;
      if ( dumi < fracSample2 ) {
	++n2;
	sow.write(b2);
      }
      ++t2;      
      b2 = sr2.read() ? sr2.cursor() : NULL;
    }

    // decide whether to sample the read or not
    humi = (str_hash(bam_name(b1)) % UINT16_MAX);
    dumi = (double)(humi+0.5) / (double)UINT16_MAX;
    if ( dumi < fracSample1 ) {
      ++n1;
      sow.write(b1);
    }
    ++t1;
  }

  notice("Closing the file");

  while ( !sr2.eof ) {
    b2 = sr2.cursor();
    humi = (str_hash(bam_name(b2)) % UINT16_MAX);
    dumi = (double)(humi+0.5) / (double)UINT16_MAX;
    if ( dumi < fracSample2 ) {
      ++n2;
      sow.write(b2);
    }
    ++t2;    
    sr2.read();
  }  

  sr1.close();
  sr2.close();  
  sow.close();

  notice("Finished writing BAM file with %d/%d (%.5lf) and %d/%d (%.5lf) reads from each file", n1, t1, (double)n1/(double)t1, n2, t2, (double)n2/(double)t2);
  */
  return 0;
}
