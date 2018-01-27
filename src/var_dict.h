#ifndef __VAR_DICT_H
#define __VAR_DICT_H

#include <cstring>
#include <vector>
#include <map>
#include <algorithm>

#include "Error.h"

struct var_dict_key {
  int64_t pos1;
  std::string ref;
  std::vector<std::string> alts;

  bool operator< (const struct var_dict_key& other) const {
    if ( pos1 < other.pos1 ) return true;
    else if ( pos1 == other.pos1 ) {
      if ( ref < other.ref ) return true;
      else if ( ref == other.ref ) {
	if ( alts.size() < other.alts.size() ) return true;
	else if ( alts.size() == other.alts.size() ) {
	  for(size_t i = 0; i < alts.size(); ++i) {
	    if ( alts[i] < other.alts[i] ) return true;
	    else if ( alts[i] > other.alts[i] ) return false;
	  }
	  return false;
	}
	else return false;
      }
      else return false;
    }
    else return false;
  }
};

typedef struct var_dict_key var_dict_key_t;

template <class T>
class var_elem {
 public:
  std::string chr;
  var_dict_key var;
  T val;
 var_elem(const char* _chr, const var_dict_key& _var, const T& _val) :
  chr(_chr), var(_var), val(_val) {}
};

typedef struct chr_var_plp chr_var_plp_t;

template <class T>
class var_dict {
 public:
  typedef typename std::map< std::string, std::map<var_dict_key_t, T> >::iterator var_dict_it_t;
  typedef typename std::map<var_dict_key_t, T>::iterator var_it_t;
  int32_t        nvar;
   
 protected:
  std::map< std::string, std::map<var_dict_key_t, T> > dict;
  std::string    tmp_chr;
  var_dict_key_t tmp_key;
  var_it_t       tmp_it;

 public:
 var_dict() : nvar(0) {}
  static bool parse_var(const char* varid, std::string& chr, var_dict_key_t& key) {
    const char* s = varid;
    
    while( ( *s ) && ( *s != ':' ) ) ++s;
    if ( *s == '\0' ) return false;
    chr.assign(varid, s-varid);
    s = varid = s+1;

    while( ( *s ) && ( *s != ':' ) ) ++s;
    if ( *s == '\0' ) return false;
    key.pos1 = (uint64_t)atoll(varid);
    s = varid = s+1;

    while( ( *s ) && ( *s != '_' ) ) ++s;
    if ( *s == '\0' ) return false;
    key.ref.assign(varid, s-varid);
    s = varid = s+1;

    for(int32_t i=0; *s; ++i) {
      while( ( *s ) && ( *s != ',' ) ) ++s;
      key.alts.resize(i+1);
      key.alts[i].assign(varid, s-varid);
      if ( *s == '\0' ) return true;
      s = varid = s+1;
    }
    return false;
  }

  T& operator[]( const char* varid ) {
    if ( !parse_var(varid, tmp_chr, tmp_key) )
      error("[E:%s:%d %s] Cannot parse variant id %s", __FILE__, __LINE__, __FUNCTION__, varid);
    if ( dict[tmp_chr].find(tmp_key) == dict[tmp_chr].end() ) {
      ++nvar;
    }
    return dict[tmp_chr][tmp_key];
  }
  
  bool has_var( const char* varid ) {
    if ( !parse_var(varid, tmp_chr, tmp_key) )
      error("[E:%s:%d %s] Cannot parse variant id %s", __FILE__, __LINE__, __FUNCTION__, varid);
    return ( (tmp_it = dict[tmp_chr].find(tmp_key)) != dict[tmp_chr].end() );
  }

  //var_dict_it_t get_dict_begin() const { tmp_chr = dict.begin()->first; return dict.begin(); }
  //var_dict_it_t get_dict_end() const { tmp_chr = ""; return dict.end(); }

  //var_it_t get_cur_it() const { return tmp_it; }

  var_it_t get_chr_beg(const char* chr) { return ( tmp_it = dict[(tmp_chr = chr)].begin() ); }
  var_it_t get_chr_end(const char* chr) { return ( tmp_it = dict[(tmp_chr = chr)].end() ); }

  //var_it_t get_next_it() { return tmp_it++; }

  var_it_t lower_bound(const char* varid ) {
    if ( !parse_var(varid, tmp_chr, tmp_key) )
      error("[E:%s:%d %s] Cannot parse variant id %s", __FILE__, __LINE__, __FUNCTION__, varid);
    return ( dict[tmp_chr].lower_bound(tmp_key) );
  }

  const char* get_cur_chrom() const { return tmp_chr.c_str(); }

  var_it_t lower_bound(const char* chr, uint64_t pos1) {
    tmp_key.pos1 = pos1+1;
    tmp_key.ref.clear();
    tmp_key.alts.clear();
    return ( dict[tmp_chr = chr].lower_bound(tmp_key) );
  }
  
  var_it_t upper_bound( const char* varid ) {
    if ( !parse_var(varid, tmp_chr, tmp_key) )
      error("[E:%s:%d %s] Cannot parse variant id %s", __FILE__, __LINE__, __FUNCTION__, varid);    
    return ( dict[tmp_chr].upper_bound(tmp_key) );    
  }

  var_it_t upper_bound(const char* chr, uint64_t pos1) {
    tmp_key.pos1 = pos1;
    tmp_key.ref.clear();
    tmp_key.alts.clear();
    return ( dict[tmp_chr = chr].upper_bound(tmp_key) );
  }

  int32_t get_depth_distribution(std::vector<int32_t>& v) {
    var_dict_it_t d_it = dict.begin();
    var_dict_it_t d_it_tmp;
    while( d_it != dict.end() ) {
      notice("Calculating depth distribution for %s", d_it->first.c_str());
      var_it_t it = d_it->second.begin();
      while( it != d_it->second.end() ) {
	int32_t d = it->second.depth();
	if ( d >= (int32_t)v.size() )
	  v.resize(d+1,0);
	++v[d];
	++it;
      }
      ++d_it;
    }
    return (int32_t)v.size();    
  }

  int32_t vectorize(std::vector< var_elem<T> >& v, bool remove = true) {
    var_dict_it_t d_it = dict.begin();
    var_dict_it_t d_it_tmp;
    while( d_it != dict.end() ) {
      notice("Processing %s", d_it->first.c_str());
      var_it_t it = d_it->second.begin();
      while( it != d_it->second.end() ) {
	v.push_back( var_elem<T>(d_it->first.c_str(), it->first, it->second) );
	++it;
      }

      if ( remove ) {
	d_it_tmp = d_it;
	++d_it;
	dict.erase(d_it_tmp);
      }
    }
    return (int32_t)v.size();
  }
};

#endif
