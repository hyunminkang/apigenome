#ifndef __TSV_READER_H
#define __TSV_READER_H

#include <cstdlib>
#include <cstring>
#include <vector>
#include "htslib/kstring.h"
#include "htslib/kseq.h"
#include "htslib/hts.h"
#include "Error.h"

class tsv_reader {
public:
  htsFile* hp;
  kstring_t str;
  int32_t lstr;
  int32_t nfields;
  int32_t* fields;
  int32_t nlines;

  
  bool open(const char* filename);
  int32_t read_line();
  const char* str_field_at(int32_t idx);
  int32_t int_field_at(int32_t idx);
  double double_field_at(int32_t idx);
  int32_t store_to_vector(std::vector<std::string>& v);

  tsv_reader() : hp(NULL), lstr(0), nfields(0), fields(NULL), nlines(0) {
    str.l = str.m = 0; str.s = NULL;
  }

  tsv_reader(const char* filename) : lstr(0), nfields(0), fields(NULL), nlines(0) {
    str.l = str.m = 0; str.s = NULL;    
    open(filename);
  }  

  ~tsv_reader() {
    if ( str.s ) free(str.s);
    if ( fields ) free(fields);
  }
};
#endif
