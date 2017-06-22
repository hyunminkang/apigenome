#ifndef __APIGENOME_FASTA_GC_H
#define __APIGENOME_FASTA_GC_H

#include <string>
#include <vector>
#include <map>
#include <cstdlib>
#include <stdint.h>

extern "C" {
#include "htslib/hts.h"
#include "htslib/faidx.h"
#include "htslib/bgzf.h"
}

// GC value will be 0-65534, with 65535 be the region with N in the targeted base
class fastaGC {
public:
  faidx_t* fai;
  BGZF* fp;
  uint16_t sliding_unit; 
  uint16_t window_size;
  int32_t cur_tid;
  int32_t cur_pos1;
  uint16_t cur_gc;
  bool on_memory;

  std::vector<std::string> seqnames;
  std::vector<int64_t> seqlens;
  std::vector<int64_t> offsets;
  std::vector<uint16_t*> mem_gcs;
  std::map<std::string,int32_t> seq2tid;

 fastaGC() : fai(NULL), fp(NULL), sliding_unit(0), window_size(0), cur_tid(-1), cur_pos1(-1), cur_gc((uint16_t)65535u), on_memory(false) {}
  ~fastaGC();
  
  bool openFasta(const char* filename);
  bool writeGC(const char* filename, uint16_t window, uint16_t slide);
  bool openGC(const char* filename);
  
  uint16_t getGC(const char* chr, int64_t pos1);
  uint16_t nextGC(int32_t pos_to_add = 1);
};

#endif 
