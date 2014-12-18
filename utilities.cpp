


#include "utilities.h"
#include "matrix.h"
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <iostream>

using namespace std;



void cal_niw_posterior(matrix pri_mu0, double pri_m, matrix pri_psi, double pri_n0, vector<matrix> lik_x, matrix& pos_mu0, double& pos_m, matrix& pos_psi, double& pos_n0)
{
    pos_m = pri_m + lik_x.size();
    pos_n0 = pri_n0 + lik_x.size();
    pos_mu0 = pri_mu0;
    pos_psi = pri_psi;
    if (lik_x.size()) {
        matrix xba = matrix(lik_x[0].m, 1);
        for (int i = 0; i < lik_x.size(); i++) {
            xba = xba + (lik_x[i] * (1.0 / lik_x.size()));
        }
        matrix S = matrix(lik_x[0].m, lik_x[0].m);
        for (int i = 0; i < lik_x.size(); i++) {
            S = S + ((lik_x[i] - xba) * (lik_x[i] - xba).transpose());
        }
        pos_mu0 = ((xba * lik_x.size()) + (pri_mu0 * pri_m)) * (1.0 / pos_m);
        pos_psi = pri_psi + S + ((xba - pri_mu0) * (xba - pri_mu0).transpose()) * (lik_x.size() * pri_m / double(pos_m));
    }
}


double cal_kl_normal(matrix mu0, matrix sigma0, matrix mu1, matrix sigma1) {
    return 0.5 * ((sigma1.inverse() * sigma0).trace() + ((mu1 - mu0).transpose() * sigma1.inverse() * (mu1 - mu0)).num[0][0] - log(sigma0.det() / sigma1.det()));
}




