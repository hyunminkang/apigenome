#ifndef __CRAMORE_H
#define __CRAMORE_H

#include <getopt.h>
#include <cstdio>
#include <vector>
#include <algorithm>
#include <utility>
#include <string>
#include <set>

// Hyun's codes
#include "params.h"
#include "Error.h"
#include "PhredHelper.h"

// Adrian's codes
#include "genome_interval.h"
#include "bcf_ordered_reader.h"
#include "bcf_ordered_writer.h"
#include "hts_utils.h"

// bcftools's code
#include "filter.h"

// htslib's code


// others
struct _bcf_vfilter_arg {
  std::vector<std::string> required_filters; // require at least one of the
  std::string include_expr;
  std::string exclude_expr;
};

struct _bcf_gfilter_arg {
  int32_t minDP;
  int32_t minGQ;
  //int32_t minAD;

  _bcf_gfilter_arg() {
    minDP = 0;
    minGQ = 0;
  }
};

typedef struct _bcf_vfilter_arg bcf_vfilter_arg;
typedef struct _bcf_gfilter_arg bcf_gfilter_arg;

#define MASK_GT_MISS   0x01
#define MASK_GT_HOMREF 0x02
#define MASK_GT_HET    0x04
#define MASK_GT_HOMALT 0x08
#define MASK_GT_NONREF (MASK_GT_HET|MASK_GT_HOMALT)
#define MASK_GT_NOMISS (MASK_GT_HOMREF|MASK_GT_HET|MASK_GT_HOMALT)
#define MASK_GT_ALL    (MASK_GT_MISS|MASK_GT_HOMREF|MASK_GT_HET|MASK_GT_HOMALT)

#define FLT_INCLUDE 1
#define FLT_EXCLUDE 2



#endif
