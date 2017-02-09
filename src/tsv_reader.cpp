#include "tsv_reader.h"

/*
  htsFile* hp;
  kstring_t str;
  int32_t lstr;
  int32_t nfields;
  int32_t* fields;
  int32_t nlines;
*/

bool tsv_reader::open(const char* filename) {
  hp = hts_open(filename, "r");
  if ( hp == NULL ) return false;
  return true;
    //error("[E:%s:%s %s] Cannot open file %s for reading", __FILE__, __LINE__, __FUNCTION__, filename);
}

int32_t tsv_reader::read_line() {
  lstr = hts_getline(hp, KS_SEP_LINE, &str);
  if ( lstr <= 0 ) {
    nfields = 0;
    return lstr;
  }
  fields = ksplit(&str, 0, &nfields);
  ++nlines;
  return nfields;
}

const char* tsv_reader::str_field_at(int32_t idx) {
  if ( idx >= nfields ) {
    error("[E:%s:%s %s] Cannot access field at %d >= %d", __FILE__, __LINE__, __FUNCTION__, idx, nfields);
    //return NULL;
  }
  return ( &str.s[fields[idx]] );
}

int32_t tsv_reader::int_field_at(int32_t idx) {
  if ( idx >= nfields )
    error("[E:%s:%s %s] Cannot access field at %d >= %d", __FILE__, __LINE__, __FUNCTION__, idx, nfields);    
  return ( atoi(&str.s[fields[idx]]) );
}

double tsv_reader::double_field_at(int32_t idx) {
  if ( idx >= nfields )
    error("[E:%s:%s %s] Cannot access field at %d >= %d", __FILE__, __LINE__, __FUNCTION__, idx, nfields);    
  return ( atof(&str.s[fields[idx]]) );
}

int32_t tsv_reader::store_to_vector(std::vector<std::string>& v) {
  v.resize(nfields);
  for(int32_t i=0; i < nfields; ++i) {
    v[i].assign(&str.s[fields[i]]);
  }
  return nfields;
}
