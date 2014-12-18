

#ifndef __TRAIN_H__
#define __TRAIN_H__

#include "matrix.h"
#include "stored_info.h"
#include <stdlib.h>
#include <math.h>
#include <vector>
using namespace std;



void initialization(stored_info *info);

void sampling_cluster_parameters(stored_info *info);

void sampling_supercluster_parameters(stored_info *info);

void sampling_data_parameters(stored_info *info);

void splitting_clusters(stored_info *info);

void train(stored_info *info);

double cross_validiction(stored_info *info);



#endif
