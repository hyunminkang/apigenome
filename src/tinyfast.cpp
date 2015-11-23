#include "params.h"
#include "Error.h"
#include "wFile.h"
#include "pFile.h"

#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <cstring>

int runGlfDepthSum(int argc, char** argv) {
  paramList pl;

  int minMapQ = 0;
  bool countN = false;
  int depthCap = 255;
  int minDepth = 0;
  int beg = 0;
  int end = 0x7fffffff;
  std::string xChr("X");
  std::string yChr("Y");  
  int xStart = 2699520;
  int xStop = 154931044;
  std::string inf("/dev/stdin");
  std::string outf("/dev/stdout");
  
  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input and Output", NULL)
    LONG_STRING_PARAM("in",&inf,"Input File")
    LONG_STRING_PARAM("out",&outf,"Output File")
    
    LONG_PARAM_GROUP("Filtering bases", NULL)
    LONG_INT_PARAM("minMQ",&minMapQ, "Minimum mapping quality thresholds")
    LONG_PARAM("countN",&countN,"Do not exclude based where the reference is N")
    LONG_INT_PARAM("beg",&beg,"Base position to begin calculating depth")
    LONG_INT_PARAM("end",&end,"Base position to end calculating depth")
    LONG_INT_PARAM("depth-cap",&depthCap,"Maximum depth to be capped to")
    LONG_INT_PARAM("min-depth",&minDepth,"Minimum depth to be included as eligible")    

    LONG_PARAM_GROUP("Sex chromosome options", NULL)
    LONG_STRING_PARAM("xChr",&xChr,"X chromosome label")
    LONG_STRING_PARAM("yChr",&yChr,"Y chromosome label")
    LONG_INT_PARAM("xStart",&xStart,"Start base position of non-PAR chromosome X")
    LONG_INT_PARAM("xStop",&xStop,"Stop base position of non-PAR chromosome X")    
  END_LONG_PARAMS();

  pl.Add(new longParams("tinyfast glf-depth-sum - Depth calculator from glf file\n\nThe usage is\n\t$(GOTCLOUD_ROOT)/bin/samtools-hybrid glfview [input.glf] | tinyfast glf-depth-sum [options]\n\n* The output is the total number of eligible bases in the region (denominator) and the sum of the depths across the region (numerator)\n* When reading sex chromosomes, it calculates depth excluding PAR", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  // take the input from STDIN
  pFile pf(inf.c_str());
  wFile wf(outf.c_str());  
  const char* line = NULL;
  const char* ws[6];
  int l, i;
  int64_t numBases = 0;
  int64_t sumDepth = 0;
  
  for( l=0; ( line = pf.getLine() ) != NULL; ++l ) {
    if ( !isdigit(line[0]) && !isalpha(line[0]) ) continue;
    
    ws[0] = line-1;
    for(i=0; i < 5; ++i) {
      ws[i+1] = strchr(ws[i]+1,'\t');
    }

    if ( l == 0 ) {  // check chromosome
      if ( ( ws[1] - ws[0] ==  (int)xChr.size() + 1 ) && ( xChr.compare(0, ws[1]-ws[0]-1, line) == 0 ) ) {
	beg = xStart;
	end = xStop;
      }
    }

    if ( ( beg > 0 ) || ( end < 0x7fffffff ) ) {  // check position;
      int pos = atoi(&ws[1][1]);
      if ( ( pos < beg ) || ( pos > end ) )
	continue; // skip the position
    }

    if ( ws[2][1] == '*' ) continue;
    if ( ( !countN ) && toupper(ws[2][1]) == 'N' ) continue;

    if ( ( minMapQ > 0 ) && ( atoi(&ws[4][1]) < minMapQ ) ) continue;

    int depth = atoi(&ws[3][1]);
    
    if ( depth < minDepth ) continue;

    if ( depth > depthCap )
      depth = depthCap;

    ++numBases;
    sumDepth += depth;
  }

  wf.printf("%lld\t%lld\t%.4lf\n",numBases,sumDepth,(double)sumDepth/(double)numBases);

  return 0;
}

int main(int argc, char** argv) {
  if ( argc < 2 ) {
    printf("tinyfast -- Fast analytic tools for small tasks \n");
    printf("Copyright (c) 2015 Hyun Min Kang\n");
    printf("\nUsage : %s [command] [options]\n\n",argv[0]);
    printf("\tType one of the following commands below to get detailed usage\n\n");
    printf("\t%s glf-depth-sum  [options] : Calculate sum of depth piped from GLF files\n",argv[0]);
    printf("\n");
  }
  else {
    std::string cmd(argv[1]);
    if ( cmd == "glf-depth-sum" ) {
      return runGlfDepthSum(argc-1,argv+1);
    }
    else {
      error("Unrecognized command %s\n",argv[0]);
    }
  }
  return 0;
}
