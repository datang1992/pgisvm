#ifndef __CLUSTER_SAMPLING___
#define __CLUSTER_SAMPLING___


#include <math.h>
#include "sampling.h"
#include <gsl/gsl_sf.h>
#include <gsl/gsl_randist.h>

using namespace std;


void sample_cluster(stored_info *info, int index);

void sample_subcluster(stored_info *info, int index, int subindex);

void sample_supercluster(stored_info *info, int index);

void calculate_split_probability(stored_into *info, int index);

void propose_split(stored_info *info, int index);


void initial_clusers(stored_info *info);

void sampling_process(stored_info *info, int rounds);


void clean_info(stored_info *info);





#endif
