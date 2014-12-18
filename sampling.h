

#ifndef ___SAMPLING_H__
#define ___SAMPLING_H__

#include "matrix.h"
#include <stdlib.h>
#include <math.h>
#include <vector>
using namespace std;


double sampling_from_uniform_0_1();

double sampling_from_beta_1 (double alpha);

double sampling_from_gaussian (double mu, double sigma_square);

matrix sampling_from_gaussian (matrix mean, matrix cov);

double sampling_from_inverse_gaussian (double mu, double lambda);

int sampling_from_multivector(vector<double>);


#endif
