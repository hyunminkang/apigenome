#include "fastaGC.h"
#include "Error.h"
#include <cmath>

fastaGC::~fastaGC() {
  if ( fai ) {
    fai_destroy(fai);
  }
  if ( fp ) {
    if ( bgzf_close(fp) < 0 ) error("bgzf_close failed");    //bgzf_close(fp);
  }
  if ( on_memory ) {
    for(int32_t i=0; i < (int32_t)mem_gcs.size(); ++i) {
      delete [] mem_gcs[i];
    }
  }
}

bool fastaGC::openFasta(const char* filename) {
  fai = fai_load(filename);
  if ( fai == NULL ) {
    error("Cannot load fasta file %s and its genome index %s.fai", filename, filename);
    return false;
  }

  int64_t byteOffset = sizeof(int32_t)+sizeof(uint16_t)+sizeof(uint16_t);
  seqnames.resize(faidx_nseq(fai));
  seqlens.resize(faidx_nseq(fai));
  offsets.clear();
  seq2tid.clear();
  for(int32_t i=0; i < (int32_t)seqnames.size(); ++i) {
    seqnames[i] = faidx_iseq(fai, i);
    seqlens[i] = faidx_seq_len(fai, faidx_iseq(fai, i));
    seq2tid[seqnames[i]] = i;
    byteOffset += sizeof(int64_t);
    byteOffset += sizeof(int32_t);
    byteOffset += (int32_t)seqnames[i].size();
  }

  for(int32_t i=0; i < (int32_t)seqnames.size(); ++i) {
    offsets.push_back(byteOffset);
    // when the length is L, the window is W, and the bin size is B
    // the total number of bins are ceiling(L/B), with some inefficiencies
    byteOffset += sizeof(uint16_t) * ceil((double)seqlens[i]/(double)sliding_unit);
  }
  offsets.push_back(byteOffset);
  return true;
}

bool fastaGC::writeGC(const char* filename, uint16_t window, uint16_t slide) {
  // open file
  if ( fai == NULL ) {
    warning("%s failed as the fasta file was not open first", __PRETTY_FUNCTION__);
    return false;
  }

  notice("Opening %s for writing",filename);
  window_size = window;
  sliding_unit = slide;
  
  fp = bgzf_open(filename, "wb");
  if ( fp == NULL ) {
    warning("%s failed as the output file cannot be created", __PRETTY_FUNCTION__);
    return false;
  }

  if ( bgzf_index_build_init(fp) < 0 )
    error("[E:%s] Error in building BGZF index", __PRETTY_FUNCTION__);  

  if ( bgzf_write(fp, &window_size, sizeof(uint16_t)) < 0 ) return false;
  if ( bgzf_write(fp, &sliding_unit, sizeof(uint16_t)) < 0 ) return false;  
  // first 4 bytes represents the number of sequences
  int32_t nseq = (int32_t)seqnames.size();
  if ( bgzf_write(fp, &nseq, sizeof(int32_t)) < 0 ) return false;

  // write the chromosome length, and length of chromosome name
  // and the chromosome name in order
  for(int32_t i=0; i < nseq; ++i) {
    if ( bgzf_write(fp, &seqlens[i], sizeof(int64_t)) < 0 ) return false;
    int32_t len = (int32_t)seqnames[i].size();
    if ( bgzf_write(fp, &len, sizeof(int32_t)) < 0 ) return false;
    if ( bgzf_write(fp, seqnames[i].c_str(), sizeof(char)*len) < 0 ) return false;
  }

  // for each chromosome, write the GC contents in 2 bytes, normalized by the bin size
  for(int32_t i=0; i < nseq; ++i) {
    notice("Processing contig %s...",seqnames[i].c_str());
    // read the entire chromosome. don't worry about the memory
    int32_t len;
    char* seq = fai_fetch(fai, seqnames[i].c_str(), &len);
    if ( len != (int32_t)seqlens[i] )
      error("Sequence length %lld does not match with the length of fetched sequences %d", seqlens[i], len);
    int32_t ats = 0, gcs = 0, ns = 0;
    int64_t beg0 = 0;
    int64_t end0 = 0, newbeg0, newend0;
    for(int64_t pos0 = 0; pos0 < seqlens[i]; pos0 += sliding_unit) {
      //if ( pos0 % 10000 == 0 ) notice("Processing %lld",pos0);
      newend0 = pos0 + sliding_unit;
      newbeg0 = newend0 < window_size ? 0 : newend0 - window_size;
      
      for(int64_t cur0 = end0; cur0 < newend0; ++cur0) {
	switch(seq[cur0]) {
	case 'A': case 'T': case 'a': case 't':
	  ++ats;
	  break;
	case 'C': case 'G': case 'c': case 'g':
	  ++gcs;
	  break;
	default:
	  ++ns;
	  break;
	}
      }

      for(int64_t cur0 = beg0; cur0 < newbeg0; ++cur0) {
	switch(seq[cur0]) {
	case 'A': case 'T': case 'a': case 't':
	  --ats;
	  break;
	case 'C': case 'G': case 'c': case 'g':
	  --gcs;
	  break;
	default:
	  --ns;
	  break;
	}
      }

      uint16_t gc = (gcs+ats == 0) ? 65535u : (uint16_t)floor((double)gcs/(double)(gcs+ats)*window_size);
      if ( bgzf_write(fp, &gc, sizeof(uint16_t)) < 0 ) return false;

      beg0 = newbeg0;
      end0 = newend0;
    }
    free(seq);
  }

  if ( bgzf_index_dump(fp, (std::string(filename)+".gzi").c_str(), NULL) < 0 )
    error("Error in writing index file %s.gzi", filename);

  if ( bgzf_close(fp) < 0 ) error("Close failed");
  fp = NULL;
  notice("Finished writing file %s",filename);
  
  return true;
}

bool fastaGC::openGC(const char* filename) {
  fp = bgzf_open(filename, "rb");
  if ( fp == NULL ) {
    warning("%s failed as the input file %s failed to open", __PRETTY_FUNCTION__,filename);
    return false;
  }
  int32_t nseq;

  if ( bgzf_read(fp, &window_size, sizeof(uint16_t)) != sizeof(uint16_t) ) {
    error("First 2 bytes cannot be read");
    return false;    
  }

  if ( bgzf_read(fp, &sliding_unit, sizeof(uint16_t)) != sizeof(uint16_t) ) {
    error("Next 2 bytes cannot be read");
    return false;    
  }  
  
  if ( bgzf_read(fp, &nseq, sizeof(int32_t)) != sizeof(int32_t) ) {
    error("Next 4 bytes cannot be read");
    return false;
  }
  seqnames.resize(nseq);
  seqlens.resize(nseq);
  seq2tid.clear();
  offsets.clear();
  int64_t byteOffset = sizeof(int32_t)+sizeof(uint16_t)+sizeof(uint16_t);  
  for(int32_t i=0; i < nseq; ++i) {
    if ( bgzf_read(fp, &seqlens[i], sizeof(int64_t)) != sizeof(int64_t) ) {
      error("Failed reading first 8 bytes tid=%d",i);      
      return false;
    }
    int32_t len = 0;
    if ( bgzf_read(fp, &len, sizeof(int32_t)) != sizeof(int32_t) ) {
      error("Failed reading second 4 bytes tid=%d",i);            
      return false;
    }
    char* chrom = new char[len+1];
    if ( bgzf_read(fp, chrom, sizeof(char)*len) != (ssize_t)(sizeof(char)*len) ) {
      error("Failed reading %s bytes tid=%d",len,i);                  
      return false;
    }
    chrom[len] = '\0';
    seqnames[i] = chrom;
    //seqnames[i].assign(chrom);
    delete[] chrom;
    seq2tid[seqnames[i]] = i;

    byteOffset += sizeof(int64_t);
    byteOffset += sizeof(int32_t);
    byteOffset += (int32_t)seqnames[i].size();    
  }
  for(int32_t i=0; i < (int32_t)seqnames.size(); ++i) {
    offsets.push_back(byteOffset);
    // when the length is L, the window is W, and the bin size is B
    // the total number of bins are ceiling(L/B), with some inefficiencies
    int32_t nbins = (int32_t)ceil((double)seqlens[i]/(double)sliding_unit);    
    byteOffset += ( sizeof(uint16_t) * nbins );
    if ( on_memory ) {
      mem_gcs.push_back( new uint16_t[nbins] );
      if ( bgzf_read(fp, mem_gcs.back(), nbins * sizeof(uint16_t) ) != (ssize_t)(sizeof(uint16_t)*nbins) ) {
	error("Failed reading stored GC contents for %s",seqnames[i].c_str());
      }
    }
  }
  offsets.push_back(byteOffset);

  if ( bgzf_index_load(fp, filename, ".gzi") < 0 )
    error("Could not load index file: %s.gzi", filename);
  
  return true;
}

uint16_t fastaGC::getGC(const char* chr, int64_t pos1) {
  if ( seq2tid.find(chr) == seq2tid.end() )
    error("[E:%s] Cannot find chromosome %s", __PRETTY_FUNCTION__, chr);
  cur_tid = seq2tid[chr];
  cur_pos1 = pos1;
  
  if ( on_memory ) {
    return (cur_gc = mem_gcs[cur_tid][(cur_pos1-1)/sliding_unit]);
  }
  else {
    int64_t offset = offsets[cur_tid] + ((cur_pos1-1)/sliding_unit)*sizeof(uint16_t);
    if ( bgzf_useek(fp, offset, SEEK_SET) < 0 )
      error("[E:%s] Cannot bgzf_useek at position %s:%lld = %lld", chr, cur_pos1, offset);
    if ( bgzf_read(fp, &cur_gc, sizeof(int16_t)) != sizeof(int16_t) )
      error("[E:%s] Cannot bgzf_read at position %s:%lld = %lld", chr, cur_pos1, offset);

    return cur_gc;
  }
}

uint16_t fastaGC::nextGC(int32_t pos_to_add) {
  if ( cur_pos1 + pos_to_add > seqlens[cur_tid] ) {
    cur_pos1 = 0;
    ++cur_tid;
    return getGC(seqnames[cur_tid].c_str(),cur_pos1);
  }
  else {
    if ( on_memory ) {
      cur_pos1 += pos_to_add;
      cur_gc = mem_gcs[cur_tid][(cur_pos1-1)/sliding_unit];
    }
    else {
      int64_t cur_offset = offsets[cur_tid] + ((cur_pos1-1)/sliding_unit)*sizeof(uint16_t);
      cur_pos1 += pos_to_add;
      int64_t new_offset = offsets[cur_tid] + ((cur_pos1-1)/sliding_unit)*sizeof(uint16_t);
      if ( new_offset == cur_offset ) return cur_gc;
      
      if ( bgzf_useek(fp, new_offset - cur_offset, SEEK_CUR) < 0 )
	error("[E:%s] Cannot bgzf_useek at position %s:%lld = %lld", seqnames[cur_tid].c_str(), cur_pos1, new_offset);
      
      if ( bgzf_read(fp, &cur_gc, sizeof(int16_t)) != sizeof(int16_t) )
	error("[E:%s] Cannot bgzf_read at position %s:%lld = %lld", seqnames[cur_tid].c_str(), cur_pos1, new_offset);  
    
    }
    return cur_gc;
  }
}
