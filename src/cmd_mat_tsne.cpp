#include "cramore.h"

/* cmd_mat_tsne.cpp */

#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

int32_t cmdMatTSNE(int32_t argc, char** argv) {
  std::string input_file_name;
  std::string output_prefix;

  int32_t ndims = 2;
  double theta = 0;
  double preplexity = 50.;
  int32_t seed = 0;
  
  paramList pl;
  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Required parameters", NULL)
    LONG_STRING_PARAM("in",&input_file_name, "Input file name (in dense or sparse matrix format)")
    LONG_STRING_PARAM("out",&output_prefix, "Output file prefix")    

    LONG_PARAM_GROUP("TSNE parameters", NULL)
    LONG_INT_PARAM("dim",&ndims, "Number of dimensions to embed to")
    LONG_DOUBLE_PARAM("theta",&theta, "Theta parameter")    
    LONG_DOUBLE_PARAM("perplexity",&perplexity, "Perplexity parameter")
    LONG_INT_PARAM("seed",&seed, "Random seed")       
   END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();
  
  // Define some variables
  int32_t origN, N, D, *landmarks;
  double perc_landmarks;
  double *data;
  
  TSNE* tsne = new TSNE();
  
  // Read the parameters and the dataset
  if(tsne->load_matrix(filename, theta, perplexity, no_dims, &data, &origN, &D, &rand_seed)) {
    // Make dummy landmarks
    N = origN;
    int* landmarks = (int*) malloc(N * sizeof(int));
    if(landmarks == NULL) { printf("Memory allocation failed!\n"); exit(1); }
    for(int n = 0; n < N; n++) landmarks[n] = n;
    
    // Now fire up the SNE implementation
    double* Y = (double*) malloc(N * no_dims * sizeof(double));
    double* costs = (double*) calloc(N, sizeof(double));
    if(Y == NULL || costs == NULL) { printf("Memory allocation failed!\n"); exit(1); }
    tsne->run(data, N, D, Y, no_dims, perplexity, theta, rand_seed, false);
    
    // Save the results
    tsne->save_matrix(filename, Y, landmarks, costs, N, no_dims, (argc < 3) ? "0" : argv[2]);
    
    // Clean up the memory
    free(data); data = NULL;
    free(Y); Y = NULL;
    free(costs); costs = NULL;
    free(landmarks); landmarks = NULL;
  }
  delete(tsne);
}
