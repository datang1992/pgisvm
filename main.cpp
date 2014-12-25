#include "sampling.h"
#include "train.h"
#include "stored_info.h"
#include "time.h"
#include <gsl/gsl_sf.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_cdf.h>
#include <stdlib.h>
#include <iostream>
using namespace std;


int main()
{
	srand(time(0));
    
    /*
	#pragma	omp parallel for
	for (int i = 0;  i< 100; i++)
	{
	//	#pragma omp critical
		{
		cout << sampling_from_gaussian(0, 1)  << endl;
		}
	}

	cout << gsl_sf_gamma(0.5) << endl;

	gsl_rng *rand_gen;
	rand_gen = gsl_rng_alloc(gsl_rng_rand48);
	long seed = clock ();
	gsl_rng_set(rand_gen, seed);
	
	double alpha[3] = {3,3,3};
	double theta[3];
	
	for(int i = 0; i < 3; i++)
	{
		gsl_ran_dirichlet(rand_gen, 3, alpha, theta);

		cout << theta[0] << ' ' << theta[1] << ' ' << theta[2] << endl;
	}

	gsl_rng_free(rand_gen);
	*/
    
    int N = 1000, dim = 2;
    matrix mu = matrix(dim, 1), sigma = matrix(dim, -2);
    
    stored_info *info = new stored_info(N, dim);

    
    for (int i = 0; i < N; i++) {
        info -> x.push_back(sampling_from_gaussian(mu, sigma));
        if (-info -> x[i].num[0][0] +  info -> x[i].num[1][0] >= 0.5)
                info -> y.push_back(1);
        else
                info -> y.push_back(-1);
        
        //cout << "(" << info -> x[i].num[0][0] << ", " << info -> x[i].num[1][0] << ")  " << info -> y[i] << endl;
    }
    
    cout << "Accuracy: " << cross_validiction(info) << endl;
    cout << "Number of clusters: " << info -> K << endl;
    
    for (int i = 0; i < info -> K; i++) {
        if (info -> Num[i]) {
            cout << "Classifier " << i << ": Num = " << info -> Num[i] << ", eta = (" << info -> eta[i].num[0][0] << ", " << info -> eta[i].num[1][0] << ")" << ", b = " << info -> b[i] << endl;
        }
    }
    
    return 0;
}


