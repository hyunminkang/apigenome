#include "cramore.h"
#include "bam_ordered_reader.h"

class tssPeak {
public:
  int32_t pos;
  double intensity;
  tssPeak(int32_t _pos, double _intensity) : pos(_pos), intensity(_intensity) {}
};

class tssPeakWindow {
public:
  int32_t window;
  int32_t curTid;
  int32_t curPos;

  std::vector<tssPeak> plus;
  std::vector<tssPeak> minus;

  tssPeakWindow(int32_t win) : window(win * 2 + 1), curTid(-1), curPos(-1) {}

  void addPeak(int32_t tid, int32_t pos, double intensity, bool fwd) {
    if ( intensity >= 3 ) {
      notice("Adding peak %d:%d %lg %d", tid, pos, intensity, fwd);
    }
    //moveCurPosTo(tid, pos);
  }
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

// procapWindow class (could be moved to separate header files later on)
// This contains 5'end sites from each strand of procap reads, up to a certain window size
// It also maintains a structure of bipartite graph 
class proCapWindow {
public:
  int32_t window;
  int32_t curTid;
  int32_t curPos;
  
  tssPeakWindow tpw;

  std::vector<proCapBase> plus;
  std::vector<proCapBase> minus;
  
  proCapWindow(int32_t peakWin, int32_t pairWin) : window(peakWin * 2 + 1), curTid(-1), curPos(-1), tpw(pairWin) {
    plus.resize(window);
    minus.resize(window);
  };

  /*
  bool moveCurPosTo(int32_t tid, int32_t pos) {
    // need to lcean from from curPos-window+1 to pos-window
    if ( curTid < tid ) {
      curTid = tid;
      for(int32_t i=0; i < window; ++i) {
	plus[i].clear();
	minus[i].clear();
      }
      curPos = pos;
      return true;
    }
    else if ( curTid == tid ) {
      if ( curPos < pos ) {
	for(int32_t i=(curPos > window ? curPos-window+1 : 0); (i <= pos-window) && (i < curPos + window); ++i) {
	  plus[i % window].clear();
	  minus[i % window].clear();
	}
	curPos = pos;
	return true;
      }
      else if ( pos > curPos - window ) { // lagging behind but still valid
	return false;
      }
    }

    return false;
  }
  */

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
    
    int32_t halfWin = (window-1)/2;
    std::map<int32_t,int32_t> histPlus, histMinus;
    for(int32_t i=curPos-window+1; i <= curPos; ++i) {
      ++histPlus[plus[i % window].uniqCount()];
      ++histMinus[minus[i % window].uniqCount()];      
    }

    int32_t endPos = (pos > curPos+window) ? curPos+window : pos;

    for(int32_t i=curPos; i < endPos; ++i) {
      if ( ( plus[(i - halfWin) % window].uniqCount() == histPlus.rbegin()->first ) && ( histPlus.rbegin()->second == 1 ) ) {
	tpw.addPeak( curTid, i-halfWin, (double)plus[ (i-halfWin) % window ].uniqCount(), true );
      }
      
      if ( ( minus[(i - halfWin) % window].uniqCount() == histMinus.rbegin()->first ) && ( histMinus.rbegin()->second == 1 ) ) {
	tpw.addPeak( curTid, i-halfWin, (double)minus[ (i-halfWin) % window ].uniqCount(), false );
      }

      int32_t pc = plus[(i+1) % window].uniqCount();
      if ( histPlus[pc] == 1 ) histPlus.erase(pc);
      else --histPlus[pc];

      int32_t mc = minus[(i+1) % window].uniqCount();
      if ( histMinus[mc] == 1 ) histMinus.erase(mc);
      else --histMinus[mc];      

      plus[(i+1) % window].clear();
      minus[(i+1) % window].clear();
    }
  }

  // add a proCap read
  // If the same position/strand has already a read, check dup for umi and add one if needed
  // If the position is empty, create an object
  bool addRead(int32_t tid, int32_t pos, bool fwd, const std::string& umi) {
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
    
    if ( fwd ) { return plus[pos % window].add(umi); }
    else { return minus[pos % window].add(umi); }
  }

  void print() {
    printf("** %d:%d", curTid, curPos);
    for(int32_t i=(curPos > window ? curPos-window+1 : 0); i <= curPos; ++i) {
      int32_t j = i % window;
      if ( plus[j].totalCount + minus[j].totalCount > 0 ) {
	printf(" [%d:%u:%u:%u:%d]", i , plus[j].uniqCount(), minus[j].uniqCount(), plus[j].totalCount, minus[j].totalCount);
      }
    }
    printf("\n");
  }
};

int32_t cmdCramProcapDetect(int32_t argc, char** argv) {
  std::string inSam;
  std::string region;
  std::string outPrefix;
  int32_t pairWindow = 1000; // +- 1kb as default window to match pairs
  int32_t peakWindow = 250;  // +- 250bp as default window to declare a peak

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Options for input SAM/BAM/CRAM", NULL)
    LONG_STRING_PARAM("sam",&inSam, "Input SAM/BAM/CRAM file. Must be sorted by coordinates and indexed")
    LONG_STRING_PARAM("region",&region, "Region to focus on")
    LONG_INT_PARAM("pair-win",&pairWindow, "Window size (on each end) to expand to find a pair of best-matching peaks")
    LONG_INT_PARAM("peak-win",&peakWindow, "Window size (on each end) to expand to declare a peak")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out",&outPrefix,"Out prefix")
  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  if ( inSam.empty() ) {
    error("Missing required option --sam");
  }

  std::vector<GenomeInterval> intervals;
  if ( !region.empty() ) {
    parse_intervals(intervals, "", region);
  }
  BAMOrderedReader odr(inSam, intervals);
  bam1_t *b = bam_init1();

  proCapWindow pcw(peakWindow, pairWindow);
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
    pcw.addRead( b->core.tid, b->core.pos, fwd, umi );
  }
  pcw.resolvePeaks(INT_MAX);

  odr.close();

  return 0;
}
