#include "Error.h"
#include "PhredHelper.h"

#include "htslib/sam.h"

class bamPool {
 public:
  samFile* fp;
  bam_hdr_t* hdr;
  std::vector<bam1_t*> readPools;
  bam1_t* head;
  int32_t cursor;

  bamPool(samFile* _fp, bam_hdr_t* _hdr) : fp(_fp), hdr(_hdr), head(NULL), cursor(-1) {
  }

  int32_t init() {
    head = bam_init1();
    bam1_t* tmp = bam_init1();
    int32_t ret = sam_read1(fp, hdr, tmp);
    readPools.push(tmp);
    cursor = 0;
    return ret;
  }

  int32_t pushPool() {
    int32_t ret = sam_read1(fp, hdr, head);
    if ( ret < 0 ) return ret;

    while( hasSameReadName(head, readPools[cursor]) ) {
      if ( cursor + 1 == (int32_t)readPools.size() ) {
	readPools.push_back(bam_init1());
      }
      bam1_t* tmp = readPools[cursor+1];
      readPools[cursor+1] = head;
      head = tmp;
      ret = am_read1(fp, hdr, head);      
    }
    return ret
  }
    
  int32_t read() {
    return sam_read1(fp, hdr, head);
  }

  bool isSamePair() {
    if ( readPools.empty() ) return false;
    
  }
};
