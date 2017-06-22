#include "fastaMask.h"
#include "Error.h"

fastaMask::fastaMask(const char* fastaFile) {
  fai = fai_load(fastaFile);
  if ( fai == NULL )
    error("Cannot load fasta file %s and its genome index %s.fai", fastaFile, fastaFile);

  // parse the index information
  seqnames.resize(faidx_fetch_nseq(fai));
  seqlens.resize(faidx_fetch_nseq(fai));
  seqs.resize(faidx_fetch_nseq(fai));  
  for(int32_t i=0; i < seqnames.size(); ++i) {
    seqnames[i] = faidx_iseq(i);
    seqlens[i] = faidx_seq_len(fai, faidx_iseq(i));
    seqs[i] = (char*)malloc(sizeof(char)*seqlens[i]);
    seq2tid[seqnames[i]] = i;
  }

  kstring_t line = {0, 0, NULL};
  int32_t l, offset;
  int nseq = 0;
  while( ( l = bgzf_getline(fai->bgzf, '\n', &line) ) > 0 ) {
    if ( line.s[0] == '>' ) {
      if ( ( nseq > 0 ) && ( offset != seqlens[nseq-1] ) )
	error("[E:%s:%d %s] Offset = %d != seqlens[%d] = %d",__FILE__,__LINE__,__FUNCTION__, offset, nseq-1, seqlens[nseq-1]);	
	   
      offset = 0;
      if ( ( strncmp(line.s, fai->name[nseq], strlen(fai->name[nseq])) == 0 ) && isspace(line.s[strlen(fai->name[nseq])]) ) 
	++nseq;
      else
	error("[E:%s:%d %s] Chromosome line %s does not match to the expected chromosome name %s",__FILE__,__LINE__,__FUNCTION__, line.s, fai->name[nseq]);
    }
    else {
      memcpy(seqs[nseq-1]+offset,line.s,l);
      offset += l;
    }
    if ( offset > seqlens[nseq-1] )
      error("[E:%s:%d %s] Offset = %d > seqlens[%d] = %d",__FILE__,__LINE__,__FUNCTION__, offset, nseq-1, seqlens[nseq-1]);
  }
}

void fastaMask::maskRegion(const char* chrom, int32_t beg1, int32_t end0, char c) {
  std::map<std::string,int32_t>::iterator it = seq2tid.find(chrom);
  if ( it == seq2tid.end() )
    error("[E:%s:%d %s] Cannot find chromosome name %s",__FILE__,__LINE__,__FUNCTION__, chrom);
  maskRegion(it->second, beg1, end0, c);
}

void fastaMask::maskRegion(int32_t tid, int32_t beg1, int32_t end0, char c) {
  for(int32_t i=beg1-1; i < end0; ++i)
    seqs[i] = c;
}

char fastaMask::getMask(const char* chrom, int32_t pos1) {
  std::map<std::string,int32_t>::iterator it = seq2tid.find(chrom);
  if ( it == seq2tid.end() )
    error("[E:%s:%d %s] Cannot find chromosome name %s",__FILE__,__LINE__,__FUNCTION__, chrom);
  return getMask(it->second, pos1);
}

char fastaMask::getMask(int32_t tid, int32_t pos1) {
  return seqs[pos1-1];
}

void fastaMask::maskVcf(const char* vcfFile) {
  std::vector<GenomeInterval> intervals;
  BCFOrderedReader odr(vcfFile, intervals);

  vector<int32_t> rid2tid;

  bcf1_t* iv = bcf_init();
  while( odr.read(iv) ) {
    int32_t rid = iv->rid;
    if ( ( rid < rid2tid.size() ) || ( rid2tid[rid] < 0 ) ) {
      rid2tid.resize(rid+1,-1);
      const char* rname = bcf_hdr_id2name(odr.hdr, rid);
      std::map<std::string,int32_t>::iterator it = seq2tid.find(rname);
      if ( it == seq2tid.end() )
	error("[E:%s:%d %s] Cannot find chromosome name %s",__FILE__,__LINE__,__FUNCTION__, chrom);
      rid2tid[rid] = it->second;
    }
    int32_t tid = rid2tid[rid];
    //char ref = iv->d.allele[0]
    maskRegion(tid, iv->pos, iv->pos+iv->rlen-1);
  }

  odr.close();
}
