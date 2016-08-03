#ifndef __HTS_FILE__H
#define __HTS_FILE__H

#include <zlib.h>
#include <string>
#include <cstdlib>
#include <cstdio>

#include "htslib/hts.h"
#include "htslib/hfile.h"
#include "Error.h"

// This is a generalized version of bcf_ordered_reader class developed by Adrian Tan
// It should be able to handle BAM, CRAM, VCF, BCF, BED files
class hFile {
 public:
  std::string file_name;
  hFILE* fp;
  htsFormat fmt;

  hts_idx_t *idx;
  hts_itr_t *itr;
  
  tbx_t *tbx;

  bcf_hdr_t *bcf_hdr;
  bam_hdr_t *bam_hdr;

  bcf1_t* bcf_record;
  bam1_t* bam_record;
  kstring_t line;
  

  bool is_binary_file;  
  bool intervals_present;
  bool index_loaded;
  bool random_access_enabled;

  std::vector<GenomeInterval> intervals;
  uint32_t interval_index;
  std::map<std::string, IntervalTree*> interval_tree;

  hFile(const char* filename, std::vector<GenomeInterval>* pIntervals = NULL, bool printHeader = false);

  bool jump_to_interval(GenomeInterval& interval);
  
  bool initialize_next_interval();

  void close();

  bool read_line();    // valid only for text format file
  bool read_bcf_record();  // valid for all types

  void load(const char* filename, const char* region = NULL, bool printHeader = false) {
    head = printHeader;
    if ( region != NULL ) reg = region;
    open(filename);
    if ( reg.empty() ) {
      head = false; // head flag is valid only with specified region
      loadAll();
    }
    else if ( head ) {  // head flag is set with region specified
      loadIndex();
      loadAll();    // load the header first, and will load the region later
    }
    else {
      loadIndex();
      loadRegion();  // load the region right away
    }
  }
  
 pFile() : head(false), iter(NULL), type(-1), line(NULL), idxconf(NULL) {
  }

  // read specifying one region
 pFile(const char* filename, const char* region = NULL, bool printHeader = false) : head(printHeader), iter(NULL), type(-1), line(NULL), idxconf(NULL) {
    load(filename, region, printHeader);
    /*
    // open the file
    if ( region != NULL ) reg = region;
    open(filename);
    if ( reg.empty() ) {
      head = false; // head flag is valid only with specified region
      loadAll();
    }
    else if ( head ) {  // head flag is set with region specified
      loadIndex();
      loadAll();    // load the header first, and will load the region later
    }
    else {
      loadIndex();
      loadRegion();  // load the region right away
    }
    */
  }

  void updateRegion(const char* region, bool sepchr = false) {
    reg = region;
    if ( reg.empty() ) {
      loadAll();
    }
    else {
      if ( idxconf == NULL ) loadIndex();
      loadRegion(sepchr);
    }
  }

  int getLength() {
    return len;
  }

  int read(void* ptr, size_t count) {
    switch(type) {
    case 0:
      return (int)fread(ptr, 1, count, fp);
    case 1: 
      return gzread(gf, ptr, count);
    case 2: // bgzipped files
      return bgzf_read(t->fp, ptr, count);
    default:
      error("pFile::read() - unknown type %d\n",type);
    }
    return 0;
  }

  const char* peekLine() { return line; }

  const char* getLine() {
    //fprintf(stderr,"gerLine() called\n");

    switch(type) {
    case 0:
      if ( line == NULL ) line = new char[MAX_LINE_SIZE];
      if ( fgets(line, MAX_LINE_SIZE, fp) != NULL ) {
	//fputs(line,stderr);
	len = strlen(line); // TODO : convert to lazy evaluation
	if ( line[len-1] == '\n' ) {
	  line[len-1] = '\0';
	  --len;
	}
      }
      else {
	if ( line != NULL ) delete [] line;
	len = 0;
	line = NULL;
      }
      return line;      
      /*
      size_t tn;
      //if ( line == NULL ) line = new char[MAX_LINE_SIZE];
      if ( (len = getline(&line, &tn, fp)) != -1 ) {
	if ( line[len-1] == '\n' ) {
	  line[len-1] = '\0';  // update carriage return to null character
	  --len;
	}
	return line;
      }
      else {
	//if ( line != NULL ) delete [] line;
	len = 0;
	return NULL;
      }
      */
    case 1:
      if ( line == NULL ) line = new char[MAX_LINE_SIZE];
      if ( gzgets(gf, line, MAX_LINE_SIZE) > 0 ) {
	len = strlen(line); // TODO : convert to lazy evaluation
	if ( line[len-1] == '\n' ) {
	  line[len-1] = '\0';
	  --len;
	}
      }
      else {
	if ( line != NULL ) delete [] line;
	len = 0;
	line = NULL;
      }
      return line;
    case 2:
      if ( iter == NULL ) return NULL;
      line = (char*)ti_read(t, iter, &len);
      if ( head ) { // if reached the end of the header
	if ( (int)(*line) != idxconf->meta_char ) { 
	  //if ( iter != NULL ) ti_iter_destroy(iter); // close existing iterator
	  head = false;
	  loadRegion();
	  return getLine();
	}
      }
      else if ( line == NULL ) {
	//fprintf(stderr,"foo\n");
	ti_iter_destroy(iter);
	iter = NULL;
	return NULL;
      }
      return line;
    default:
      warning("Attempt to read empty file unknown file type.\n");
      return NULL;
    }
  }

  bool open(const char* filename, bool forcegz = false) {  // return true if gzipped, false otherwise
    fname = filename;
    type = fileType(filename);
    if ( forcegz ) type = 1;
    switch(type) {
    case 0:
      t = NULL;
      gf = NULL;
      if ( strcmp(filename,"-") == 0 ) {
	fp = stdin;
      }
      else {
	fp = fopen(filename,"r");
      }
      return (fp != NULL);
    case 1:
      t = NULL;
      gf = gzopen(filename,"rb");
      fp = NULL;
      return (gf != NULL);
    case 2:
      //notice("open() is called");
      if ( (t = ti_open(filename,0)) == 0 ) {
	warning("Cannot open %s with tabix..\n",filename);
	return false;
      }
      gf = NULL;
      fp = NULL;
      //notice("open() is successful");
      return true;
    default:
      warning("Cannot open %s. File is not accessible\n",filename);
      return false;
      break;
    }
    if ( !reg.empty() && type < 2 ) {
      error("File %s is not indexed, so cannot be acessed with specified region",filename);
    }
  }

  void close() {
    switch( type ) {
    case 0:
      if ( fp != NULL ) fclose(fp);
      fp = NULL;
      break;
    case 1:
      if ( gf != NULL ) gzclose(gf);
      if ( line != NULL ) delete [] line;
      gf = NULL;
      line = NULL;
      break;
    case 2:
      //notice("close() is called %d %d",iter,t->idx);
      //if ( iter != NULL ) ti_iter_destroy(iter);
      if ( t != NULL ) {
	//ti_index_destroy(t->idx);
	ti_close(t);
      }
      idxconf = NULL; 
      t = NULL;
      iter = NULL;
      break;
    }
  }

  void loadAll() {
    if ( type == 2 ) {
      iter = ti_query(t,0,0,0);
    }
  }

  void loadIndex() {
    if (ti_lazy_index_load(t) < 0 ) {
      error("Failed to load the index file");
    }
    idxconf = ti_get_conf(t->idx);  
  }

  void loadRegion(bool sepchr = false) {
    //notice("ti_parse_region( %s )",reg.c_str());
    if ( iter != NULL ) ti_iter_destroy(iter); // close existing iterator
    if ( ti_parse_region(t->idx, reg.c_str(), &tid, &beg, &end) != 0 ) {
      if ( sepchr ) {
	// changes all "chrAA." to "chrBB." from the files
	std::string newfname;
	int pos = 0;
	size_t ichr = 0;
	while ( (ichr = fname.find("chr",pos)) != std::string::npos ) {
	  size_t idot = fname.find_first_of("-_./",ichr);
	  std::string newchr = reg.substr(0,reg.find(':'));
	  if ( idot == std::string::npos ) 
	    error("Cannot find '.','_','-', or '/' after chr in the filename with --sepchr option");
	  newfname += (fname.substr(pos,ichr-pos) + "chr" + newchr);
	  pos = idot;
	}
	newfname += fname.substr(pos);
	fname = newfname;

	notice("Changing the VCF file name to %s",fname.c_str());

	/*
	//notice("loadRegion(true) %s",reg.c_str());
	// assume that current filename is [prefix]chr[chr].[suffix]
	int ichr = fname.find("chr");
	int idot = fname.find('.',ichr);
	std::string newchr = reg.substr(0,reg.find(':'));
	std::string prefix = fname.substr(0,ichr);
	std::string suffix = fname.substr(idot);
	fname = prefix + "chr" + newchr + suffix;
	//notice("open(%s)",fname.c_str());
	*/

	if ( fileType(fname.c_str()) < 0 ) {
	  warning("Cannot parse region %s.. Returning empty",reg.c_str());
	  iter = NULL;
	}
	else {
	  close();
	  open(fname.c_str());
	  loadIndex();
	  loadRegion();
	}
      }
      else {
	warning("Cannot parse region %s.. Returning empty",reg.c_str());
	iter = NULL;
      }
    }
    else {
      //notice("ti_query(%x, %d, %d, %d)",t,tid,beg,end);
      iter = ti_queryi(t,tid,beg,end);
    }
  }
};

#endif // __TABIXED_FILE
