#include "cramore.h"
#include "bam_ordered_reader.h"
#include "hts_utils.h"

class procapPeak {
public:
  int32_t pos;
  int32_t baseDepth;
  double  integralDepth;
  procapPeak(int32_t _pos, int32_t _baseDepth, double _integralDepth) : pos(_pos), baseDepth(_baseDepth), integralDepth(_integralDepth) {} 
};

class proCapBase {
public:
  int32_t totalCount;
  std::map<std::string,int32_t> umiCounts;

  proCapBase() : totalCount(0) {}
  
  inline int32_t uniqCount() {
    return umiCounts.size();
  }
  
  inline bool add(const std::string& umi) {
    ++totalCount; return (++umiCounts[umi] == 1 );
  }
  
  inline void clear() {
    if ( totalCount > 0 ) {
      umiCounts.clear(); totalCount = 0; }
  }
};

class proCapQueue {
public:
  int32_t winPrev;  
  int32_t winNext;
  int32_t curTid;
  int32_t curPos;
  std::vector<proCapBase> plQueue;
  std::vector<proCapBase> mnQueue;
  
  htsFile* fpBase;
  bam_hdr_t* hdr;

  proCapQueue(int32_t _winPrev, int32_t _winNext, htsFile* _fpBase, bam_hdr_t* _hdr) {
    init(_winPrev, _winNext, _fpBase, _hdr);
  }
  
  void init(int32_t _winPrev, int32_t _winNext, htsFile* _fpBase, bam_hdr_t* _hdr) {
    curTid = -1;
    curPos = 0;
    winPrev = _winPrev;
    winNext = _winNext;
    fpBase = _fpBase;
    hdr = _hdr;

    plQueue.resize(winPrev + winNext);
    mnQueue.resize(winPrev + winNext);
  }

  void advancePos(int32_t pos) {
    int32_t endPos = ( pos > curPos + winPrev + winNext ) ? curPos + winPrev + winNext : pos;
    for(int32_t i=curPos; i < endPos; ++i) {
      int32_t b = i % (winPrev + winNext);
      
      if ( plQueue[b].uniqCount() > 0 ) {
	hprintf(fpBase, "%s\t%d\t+\t%d\t%d\n", hdr->target_name[curTid], i+1, plQueue[b].uniqCount(), plQueue[b].totalCount);
	plQueue[b].clear();
      }
      
      if ( mnQueue[b].uniqCount() > 0 ) {
	hprintf(fpBase, "%s\t%d\t-\t%d\t%d\n", hdr->target_name[curTid], i+1, mnQueue[b].uniqCount(), mnQueue[b].totalCount);
	mnQueue[b].clear();
      }      
    }
  }
  
  bool addRead(int32_t tid, int32_t pos, int32_t end5, bool fwd, const std::string& umi ) {
    //notice("addRead(%d, %d, %d, %s)",tid,pos,fwd,umi.c_str());
    if ( curTid < tid ) {
      if ( curTid >= 0 ) advancePos( INT_MAX );
      curTid = tid;
      curPos = pos;
    }
    else if ( ( curTid == tid ) && ( curPos < pos ) ) {
      advancePos( pos );
      curPos = pos;
    }
    else if ( curPos != pos ) {
      error("Cannot move back from %d:%d to %d:%d", curTid, curPos, tid, curPos);
    }

    if ( end5 > pos + winNext )
      error("The read length %d is greater than winNext=%d", end5-pos, winNext);

    if ( fwd )
      plQueue[ end5 % (winPrev + winNext) ].add(umi);
    else
      mnQueue[ end5 % (winPrev + winNext) ].add(umi);
    return true;
  }
};

/*
// procapWindow class (could be moved to separate header files later on)
// This contains 5'end sites from each strand of procap reads, up to a certain window size
// It also maintains a structure of bipartite graph 
class proCapWindow {
public:
  int32_t pairWindow;
  int32_t baseWindow;
  int32_t curTid;
  int32_t curPos;
  
  std::vector<proCapBase> plBases; // this is a window of size baseWindow where baseWindow is an odd number
  std::vector<proCapBase> mnBases;

  std::vector<procapPeak> plPeaks; // this contains all "unresolved" peaks
  std::vector<procapPeak> mnPeaks;

  bam_hdr_t* hdr;
  htsFile* fbases;
  htsFile* fpeaks;
  htsFile* fpairs;

  int32_t thresBaseDepth;
  double thresIntegralDepth;  
  
  proCapWindow(int32_t baseWidth, int32_t pairWidth, bam_hdr_t* bamHdr, htsFile* fp_bases, htsFile* fp_peaks, htsFile* fp_pairs, int32_t thresB, double thresI) : pairWindow(pairWidth * 2 + 1), baseWindow(baseWidth * 2 + 1), curTid(-1), curPos(-1), hdr(bamHdr), fbases(fp_bases), fpeaks(fp_peaks), fpairs(fp_pairs), thresBaseDepth(thresB), thresIntegralDepth(thresI) {
    plBases.resize(baseWindow);
    mnBases.resize(baseWindow);
  };

  static bool findMaxPeak(const std::vector<procapPeak>& peaks, int32_t& imax, double& max, int32_t beg, int32_t end) {
    max = -1.;
    imax = -1;
    for(int32_t i=peaks.size()-1; ( i >= 0 ) && ( peaks[i].pos >= beg ); --i) {
      if ( ( peaks[i].pos <= end) && ( max <= peaks[i].integralDepth ) ) {
	max = peaks[i].integralDepth;
	imax = i;
      }
    }
    return ( max >= 0 );
  };

  void outputPeakPair( const procapPeak& pl, const procapPeak& mn ) {
    if ( pl.pos > mn.pos ) {
      hprintf(fpairs, "%s\t%d\t%d\tdiv\t%d\t%.1lf\t%.1lf\t%d\t%d\n", hdr->target_name[curTid], mn.pos, pl.pos, pl.pos-mn.pos, mn.integralDepth, pl.integralDepth, mn.baseDepth, pl.baseDepth);
    }
    else {
      hprintf(fpairs, "%s\t%d\t%d\tconv\t%d\t%.1lf\t%.1lf\t%d\t%d\n", hdr->target_name[curTid], pl.pos, mn.pos, pl.pos-mn.pos, pl.integralDepth, mn.integralDepth, pl.baseDepth, mn.baseDepth);      
    }
  }

  // identify the best pair of peaks occurring before [pos-(window-2)/1]
  bool findBestPeakPair(int32_t pos) {
    int32_t ipl = -1, imn = -1;
    double  maxpl = -1, maxmn = -1;

    // if either peak vector is empty return no pair
    if ( !findMaxPeak(plPeaks, ipl, maxpl, 0, pos) ) return false;
    if ( !findMaxPeak(mnPeaks, imn, maxmn, 0, pos) ) return false;

    // if max is not close enough yet, return no pair
    if ( plPeaks[ipl].pos > pos - (pairWindow-1)/2 ) return false;
    if ( mnPeaks[imn].pos > pos - (pairWindow-1)/2 ) return false;

    // if max pair is close enough, select it
    if ( abs(plPeaks[ipl].pos - mnPeaks[imn].pos) < pairWindow ) {
      outputPeakPair( plPeaks[ipl], mnPeaks[imn] );
      plPeaks.erase( plPeaks.begin() + ipl );
      mnPeaks.erase( mnPeaks.begin() + imn );
      return true;
    }
    else {
      int32_t ipl2 = -1, imn2 = -1;
      double maxpl2 = -1, maxmn2 = -1;
      findMaxPeak(plPeaks, ipl2, maxpl2, mnPeaks[imn].pos-(pairWindow-1)/2, mnPeaks[imn].pos+(pairWindow-1)/2);
      findMaxPeak(mnPeaks, imn2, maxmn2, plPeaks[ipl].pos-(pairWindow-1)/2, plPeaks[ipl].pos+(pairWindow-1)/2);

      if ( ipl2 < 0 ) // minus peak does not have any valid pair, so have to be removed
	mnPeaks.erase( mnPeaks.begin() + imn );

      if ( imn2 < 0 )
	plPeaks.erase( plPeaks.begin() + ipl );

      if ( ( ipl2 < 0 ) || ( imn2 < 0 ) ) return true;
      
      else if ( plPeaks[ipl].integralDepth > mnPeaks[imn].integralDepth ) { // plus has stronger signals
	if ( abs(plPeaks[ipl].pos-mnPeaks[imn2].pos) > pairWindow )
	  error("[E:%s:%d %s] ipl-imn2 > pairWindow : ipl = %d/%d, imn = %d/%d, ipl2 = %d/%d, imn2 = %d/%d, pos=%d",
		__FILE__, __LINE__, __FUNCTION__, 
		ipl, plPeaks[ipl].pos,
		imn, mnPeaks[imn].pos,
		ipl2, plPeaks[ipl2].pos,
		imn2, mnPeaks[imn2].pos, pos );
	outputPeakPair( plPeaks[ipl], mnPeaks[imn2] );
	plPeaks.erase( plPeaks.begin() + ipl );
	mnPeaks.erase( mnPeaks.begin() + imn2 );
	return true;
      }
      else {
	if ( abs(plPeaks[ipl2].pos-mnPeaks[imn].pos) > pairWindow )
	  error("[E:%s:%d %s] ipl2-imn > pairWindow : ipl = %d/%d, imn = %d/%d, ipl2 = %d/%d, imn2 = %d/%d, pos=%d, plPeaks.size() = %u, mnPeaks.size() = %u", __FILE__, __LINE__, __FUNCTION__,
		ipl, plPeaks[ipl].pos,
		imn, mnPeaks[imn].pos,
		ipl2, plPeaks[ipl2].pos,
		imn2, mnPeaks[imn2].pos,
		pos,
		plPeaks.size(),
		mnPeaks.size() );	
	outputPeakPair( plPeaks[ipl2], mnPeaks[imn] );
	plPeaks.erase( plPeaks.begin() + ipl2 );
	mnPeaks.erase( mnPeaks.begin() + imn );
	return true;
      }
    }
  }

  void resolvePairs(int32_t pos) {
    // nothing further to insert before pos
    // before [pos-(window-2)/1], nothing is going to be paired additionally
    while( findBestPeakPair(pos) ) {
      //notice("%d",pos);
    }
  }

  void addPeak(int32_t tid, int32_t pos, int32_t baseDepth, double integralDepth, int32_t integralWidth, bool fwd) {
    if ( ( integralDepth >= thresIntegralDepth ) && ( baseDepth >= thresBaseDepth ) ) {
      //notice("Adding peak %d:%d\t%d\t%.1lf\t%s", tid, pos, baseDepth, integralDepth, fwd ? "+" : "-");
      hprintf(fpeaks, "%s\t%d\t%d\t%d\t%d\t%.1lf\t%s\n", hdr->target_name[tid], pos-(integralWidth-1)/2, pos+(integralWidth-1)/2, pos, baseDepth, integralDepth, fwd ? "+" : "-");
      //notice("Adding peak %d:%d\t%d\t%.1lf\t%s", tid, pos, baseDepth, integralDepth,

      if ( curTid < tid ) {
	resolvePairs(INT_MAX);
	curTid = tid;
	notice("Processing chromosome %s", hdr->target_name[tid]);
      }

      if ( fwd ) 
	plPeaks.push_back( procapPeak(pos, baseDepth, integralDepth ) );
      else
	mnPeaks.push_back( procapPeak(pos, baseDepth, integralDepth ) );

      resolvePairs(pos);
    }
  }

  double integrateUniqDepth(int32_t pos, int32_t width, std::vector<proCapBase>& bases) {
    int32_t peak = bases[pos % baseWindow].uniqCount();
    double sum = 0;
    for(int32_t i=1; i < width; ++i) {
      int32_t p = bases[(pos + i) % baseWindow].uniqCount();
      int32_t m = bases[(pos-i+baseWindow) % baseWindow].uniqCount();
      if ( p >= peak ) {
	//notice("foo %d %d %d",i,p,peak);    	
	return -1;
      }
      else sum += ((double)(width - i) * (double)( p + m ) / (double)width);
    }
    return sum + peak;
  }

  void resolvePeaks( int32_t pos ) {
    //notice("resolvePeaks(%d)",pos);
    // resolve peaks between [curPos-window, pos-1]
    // need to know whether (1) the value is maxima within the window
    //  this is possible when [bp-(w-1)/2, bp+(w-1)/2] is included within a region
    // after knowing which one(s) is going to be a peak or not
    // for example:
    //   curPos = 1,101
    //   pos    = 1,105
    //   win    =   101
    //   stored: 1,001 ... 1,101
    //   peakStart = 1,051
    //   peakEnd   = 1,055
    //   init : 1,001..1,101  [1,051]
    //   scan : 1,102..1,104  [1,052-1,054]
    //print();	      
    
    int32_t halfWin = (baseWindow-1)/2;

    std::map<int32_t, int32_t> plHist, mnHist;
    
    for(int32_t i=curPos-baseWindow+1; i <= curPos; ++i) {
      if ( plBases[i % baseWindow].uniqCount() > 0 )       
	++plHist[plBases[i % baseWindow].uniqCount()];

      if ( mnBases[i % baseWindow].uniqCount() > 0 )      
	++mnHist[mnBases[i % baseWindow].uniqCount()];
    }

    int32_t endPos = (pos > curPos+baseWindow) ? curPos+baseWindow : pos;

    for(int32_t i=curPos; i < endPos; ++i) {
      if ( ( !plHist.empty() ) && ( plBases[(i - halfWin) % baseWindow].uniqCount() == plHist.rbegin()->first ) ) {
	addPeak( curTid, i-halfWin, plBases[ (i-halfWin) % baseWindow ].uniqCount(), integrateUniqDepth(i-halfWin, halfWin, plBases), baseWindow, true );
      }
      
      if ( ( !mnHist.empty() ) && ( mnBases[(i - halfWin) % baseWindow].uniqCount() == mnHist.rbegin()->first ) ) {
	addPeak( curTid, i-halfWin, mnBases[ (i-halfWin) % baseWindow ].uniqCount(), integrateUniqDepth(i-halfWin, halfWin, mnBases), baseWindow, false );
      }

      int32_t pc = plBases[(i+1) % baseWindow].uniqCount();
      if ( pc > 0 ) {
	hprintf(fbases, "%s\t%d\t%d\t+\n", hdr->target_name[curTid], i+1, pc);
	if ( plHist[pc] == 1 ) plHist.erase(pc);
	else --plHist[pc];
      }

      int32_t mc = mnBases[(i+1) % baseWindow].uniqCount();
      if ( mc > 0 ) {
	hprintf(fbases, "%s\t%d\t%d\t-\n", hdr->target_name[curTid], i+1, mc);	
	if ( mnHist[mc] == 1 ) mnHist.erase(mc);
	else --mnHist[mc];
      }

      plBases[(i+1) % baseWindow].clear();
      mnBases[(i+1) % baseWindow].clear();
    }
  }

  // add a proCap read
  // If the same position/strand has already a read, check dup for umi and add one if needed
  // If the position is empty, create an object
  bool addRead(int32_t tid, int32_t pos, bool fwd, const std::string& umi ) {
    //notice("addRead(%d, %d, %d, %s)",tid,pos,fwd,umi.c_str());
    if ( curTid < tid ) {
      if ( curTid >= 0 ) resolvePeaks( INT_MAX );
      curTid = tid;
      curPos = pos;
    }
    else if ( ( curTid == tid ) && ( curPos < pos ) ) {
      resolvePeaks( pos );
      curPos = pos;
    }
    else if ( curPos != pos ) {
      error("Cannot move back from %d:%d to %d:%d", curTid, curPos, tid, curPos);
    }
    
    if ( fwd ) { return plBases[pos % baseWindow].add(umi); }
    else { return mnBases[pos % baseWindow].add(umi); }
  }

  void print() {
    printf("** %d:%d", curTid, curPos);
    for(int32_t i=(curPos > baseWindow ? curPos-baseWindow+1 : 0); i <= curPos; ++i) {
      int32_t j = i % baseWindow;
      if ( plBases[j].totalCount + mnBases[j].totalCount > 0 ) {
	printf(" [%d:%u:%u:%u:%d]", i , plBases[j].uniqCount(), mnBases[j].uniqCount(), plBases[j].totalCount, mnBases[j].totalCount);
      }
    }
    printf("\n");
  }
};
*/

int32_t cmdCramProcapDetect(int32_t argc, char** argv) {
  std::string inSam;
  std::string region;
  std::string outPrefix;
  int32_t pairWindow = 500; // +- 1kb as default window to match pairs
  int32_t peakWindow = 100;  // +- 250bp as default window to declare a peak
  int32_t thresBaseDepth = 1;
  double thresIntegralDepth = 5;

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Options for input SAM/BAM/CRAM", NULL)
    LONG_STRING_PARAM("sam",&inSam, "Input SAM/BAM/CRAM file. Must be sorted by coordinates and indexed")
    LONG_STRING_PARAM("region",&region, "Region to focus on")
    LONG_INT_PARAM("pair-win",&pairWindow, "Window size (on each end) to expand to find a pair of best-matching peaks")
    LONG_INT_PARAM("peak-win",&peakWindow, "Window size (on each end) to expand to declare a peak")
    LONG_INT_PARAM("thres-base-depth",&thresBaseDepth, "Threshold for minimum base depth to be declared as a peak")
    LONG_DOUBLE_PARAM("thres-integral-depth",&thresIntegralDepth, "Threshold for minimum integral depth to be declared as a peak")    

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out",&outPrefix,"Out prefix")
  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  if ( inSam.empty() ) {
    error("Missing required option --sam");
  }

  if ( outPrefix.empty() ) {
    error("Missing required option --out");    
  }

  htsFile* fp_base = hts_open( (outPrefix + ".bases.txt.gz").c_str(), "wz" );
  if ( fp_base == NULL )
    error("Cannot open output file %s for writing\n",(outPrefix+".bases.txt.gz").c_str());  

  /*
  htsFile* fp_peak = hts_open( (outPrefix + ".peaks.bed.gz").c_str(), "wz" );
  if ( fp_peak == NULL )
    error("Cannot open output file %s for writing\n",(outPrefix+".peaks.bed.gz").c_str());

  htsFile* fp_pair = hts_open( (outPrefix + ".pairs.bed.gz").c_str(), "wz" );
  if ( fp_peak == NULL )
    error("Cannot open output file %s for writing\n",(outPrefix+".paris.bed.gz").c_str());  
  */

  std::vector<GenomeInterval> intervals;
  if ( !region.empty() ) {
    parse_intervals(intervals, "", region);
  }
  BAMOrderedReader odr(inSam, intervals);
  bam1_t *b = bam_init1();

  proCapQueue pcq(peakWindow, peakWindow, fp_base, odr.hdr);
  //proCapWindow pcw(peakWindow, pairWindow, odr.hdr, fp_base, fp_peak, fp_pair, thresBaseDepth, thresIntegralDepth);
  std::string umi;

  while( odr.read(b) ) {
    if ( b->core.flag & 0xf04 ) continue;
    bool fwd = ( b->core.flag & 0x010 ? false : true );

    // extract barcode
    const char *prn = bam_get_qname(b);
    const char *ptmp = NULL;
    while( ( ptmp = strchr(prn, ':') ) != NULL ) {
      prn = ptmp+1;
    }
    umi.assign(prn);
    if ( fwd ) {
      pcq.addRead( b->core.tid, b->core.pos, b->core.pos, fwd, umi );
    }
    else {
      pcq.addRead( b->core.tid, b->core.pos, bam_get_end_pos1(b), fwd, umi );      
    }
  }
  pcq.advancePos(INT_MAX);

  hts_close(fp_base);
  //hts_close(fp_peak);
  //hts_close(fp_pair);  
  odr.close();  

  //if ( tbx_index_build((outPrefix + ".peaks.bed.gz").c_str(), 0, &tbx_conf_bed) )
  //  error("Failed building tabix index for %s",(outPrefix + ".peaks.bed.gz").c_str());

  return 0;
}
