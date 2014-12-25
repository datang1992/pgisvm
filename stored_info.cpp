

#include "stored_info.h"
#include "matrix.h"
#include "time.h"
#include <math.h>
#include <vector>

using namespace std;


stored_info::stored_info()
{
    
}

stored_info::stored_info(int n, int Dim)
{
	number_of_iterations = 100;
    	initial_clusters_number = 10 * log(n);
	initial_m = 1;
	initial_n0 = 1;
    	N = n;
    	dim = Dim;
    	matrix new_vec = matrix(Dim, 1);
	matrix new_ide = matrix(Dim, -2);
	initial_mu0 = new_vec;
	initial_psi = new_ide;
	C = 1;
	l = 1;
	nu = 1;
	nu2 = 1;
    	alpha = 1;
    	rng = gsl_rng_alloc(gsl_rng_rand48);
	long seed = clock();
	gsl_rng_set(rng, seed);
}

stored_info::~stored_info()
{
    gsl_rng_free(rng);
}
