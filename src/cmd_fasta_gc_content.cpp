#include "cramore.h"
#include "fastaGC.h"
#include "tsv_reader.h"

int32_t cmdFastaGCContent(int32_t argc, char** argv) {
  std::string refFasta;
  std::string gcFile;
  std::string posList;
  std::vector<std::string> positions;
  int32_t window = 150;
  int32_t sliding = 5;
  bool create = false;
  bool lookup = false;

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Running modes",NULL)
    LONG_PARAM("create",&create,"Create GC content file. Requires --ref, --prefix, --win, --slide parameters")
    LONG_PARAM("lookup",&lookup,"Lookup GC contents from a list. Requires --prefix, --list options")

    LONG_PARAM_GROUP("Input/Output files", NULL)
    LONG_STRING_PARAM("ref",&refFasta, "Indexed FASTA files to create GC content profiles from")
    LONG_STRING_PARAM("gc",&gcFile, "Name of GC content file")
    LONG_STRING_PARAM("list",&posList, "A VCF-like file that contains genomic coordinates")
    LONG_MULTI_STRING_PARAM("pos",&positions, "Genomic coordinate(s) in [CHR]:[POS] format")    

    LONG_PARAM_GROUP("Input options", NULL)
    LONG_INT_PARAM("win",&window,"Window size")
    LONG_INT_PARAM("slide",&sliding,"Unit of slidging window")
  END_LONG_PARAMS();
  
  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  if ( !create && !lookup )
    error("Either --create or --lookup option is required");
  
  // sanity check of input arguments
  if ( create ) {
    if ( refFasta.empty() || gcFile.empty() ) {
      error("with --create option, --ref and --gc parameters are required");
    }

    notice("Creating GC profile to %s with window size %d and sliding window by %d bp from %s", refFasta.c_str(), window, sliding, gcFile.c_str());
    
    fastaGC fGC;
    if ( !fGC.openFasta(refFasta.c_str()) )
      error("Cannot open FASTA file %s. Make sure that the FASTA file has valid index (.fai) file", refFasta.c_str());

    if ( !fGC.writeGC(gcFile.c_str(), (uint16_t)window, (uint16_t)sliding) )
      error("Failed creating GC content file %s. Check if you have a permission to create the files", gcFile.c_str());
  }

  if ( lookup ) {
    if ( gcFile.empty() || ( (int32_t)posList.empty() + (int32_t)positions.empty() != 1 ) ) {
      error("with --lookup option, --gc parameter is required, and either --list or --pos options are required (but not both)");
    }

    fastaGC fGC;
    if ( !fGC.openGC(gcFile.c_str()) )
      error("Failed loading GC content file %s", gcFile.c_str());

    std::vector<std::string> chrs;
    std::vector<int64_t> pos1s;
    
    if ( !posList.empty() ) {
      tsv_reader tr(posList.c_str());
      while( tr.read_line() ) {
	if ( tr.nfields == 0 ) continue;
	
	// check whether the first field has [CHR]:[POS] format
	const char* chrom = tr.str_field_at(0);
	if ( chrom[0] == '#' ) continue; // skip headers;

	if ( tr.nfields == 1 ) {
	  if ( ( strchr(chrom,':') != NULL ) ) {
	    // assume that the format is [CHR]:[POS];
	    const char* pcolon = strchr(chrom,':');
	    std::string s(chrom, pcolon-chrom);
	    chrs.push_back(s);
	    pos1s.push_back(atoll(pcolon+1));
	  }
	  else {
	    error("Cannot recognize the line %s",chrom);
	  }
	}
	else {
	  chrs.push_back(chrom);
	  pos1s.push_back(atoll(tr.str_field_at(1)));
	}
      }
    }
    else {
      for(int32_t i=0; i < (int32_t)positions.size(); ++i) {
	size_t icolon = positions[i].find(':');
	if ( icolon == std::string::npos ) {
	  error("Cannot recognize the genomic position %s",positions[i].c_str());
	}
	chrs.push_back(positions[i].substr(0,icolon));
	pos1s.push_back(atoll(positions[i].substr(icolon+1).c_str()));
      }
    }

    for(int32_t i=0; i < (int32_t)chrs.size(); ++i) {
      printf("%s\t%ld\t%u\n", chrs[i].c_str(), pos1s[i], (uint32_t)fGC.getGC(chrs[i].c_str(),pos1s[i]));
    }
  }

  return 0;
}
