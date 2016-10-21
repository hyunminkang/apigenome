#include <getopt.h>
#include <cstdio>
#include <vector>

#include "params.h"
#include "Error.h"
#include "PhredHelper.h"
#include "htslib/sam.h"
#include "genome_interval.h"
#include "hts_util.h"
#include "bcf_ordered_reader.h"
#include "filter.h"
#include "utils.h"

int32_t runVCF2GRM(int32_t argc, char** argv) {
  std::string inVcf;
  std::string 
}
