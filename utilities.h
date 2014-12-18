

#ifndef _UTILITIES_H__
#define _UTILITIES_H__


#include "matrix.h"
#include <stdlib.h>
#include <math.h>
#include <vector>

using namespace std;



void cal_niw_posterior(matrix pri_mu0, double pri_m, matrix pri_psi, double pri_n0, vector<matrix> lik_x, matrix& pos_mu0, double& pos_m, matrix& pos_psi, double& pos_n0);

double cal_kl_normal(matrix mu0, matrix sigma0, matrix mu1, matrix sigma1);




#endif
