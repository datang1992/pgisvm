


#include "sampling.h"
#include "matrix.h"
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <vector>
#include <time.h>


using namespace std;


const double pi = 3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798;


/* Using rand() function to do the sampling, four rounds.
 */
double sampling_from_uniform_0_1()
{
	double x1 = rand() % RAND_MAX, x2 = rand() % RAND_MAX, x3 = rand() % RAND_MAX, x4 = rand() % RAND_MAX;
	double x_max = RAND_MAX, x_max2 = x_max * x_max * x_max * x_max;
	return (x1 + x2 * x_max + x3 * x_max * x_max + x4 * x_max * x_max * x_max) / x_max2; 
}


/*CDF: F(x|\alpha) = 1 - (1 - x) ^ \alpha
 * Using cdf_inverse_uniform method to do the sampling/Users/tangda/Desktop/pgisvm/sampling.cpp
 */
double sampling_from_beta_1 (double alpha)
{
	double a = sampling_from_uniform_0_1();
	if(alpha <= 1e-8)
	{
		cout << "Error: Doing sampling from Beta(1,\\alpha) distribution, where \\alpha = " << alpha << " !" << endl;
		return 0.5;
	}
	double ret = 1 - pow(1 - a, 1 / alpha);
	return ret;
}

/* Sampling variables from one-dimension Gaussian Distribution.
 */
double sampling_from_gaussian (double mu, double sigma_square)
{
	if(sigma_square < 1e-8)
	{
		cout << "Error: Doing sampling from Gaussian(\\mu, \\sigma^2) distribution, where \\sigma^2 = " << sigma_square << " !" << endl;
		return mu;
	}
	double u = sampling_from_uniform_0_1(), v = sampling_from_uniform_0_1();
	double x = sqrt(-2 * log(u)) * cos(2 * pi * v);
       	return sqrt(sigma_square) * x + mu;
}

/* Sampling random vectors from a Multivariate Gaussian Distribution. Using Cholesky Decomposition.
 */
matrix sampling_from_gaussian (matrix mean, matrix cov)
{
	if(mean.n != 1 || mean.m != cov.m || mean.m != cov.n)
	{
		cout << "Error: Try to draw random variables from a Multivariate Gaussian Distribution with invalid parameters!" << endl;
		return mean;
	}

	matrix L = cov.cholesky();
	matrix ret = matrix(mean.m, 1);
	
	for(int i = 0; i < mean.m; i++)
		ret.num[i][0] = sampling_from_gaussian(0, 1);
	
	ret = mean + L * ret;
	return ret;
}

/* Sampling variables from a Inverse Gaussian Distribution.
 */
double sampling_from_inverse_gaussian (double mu, double lambda)
{
	if(mu <= 0 || lambda <= 0)
	{
		cout << "Error: Try to draw random variables from a Inverse Gaussian Distribution with parameter: \\mu = " << mu << " and \\lambda = " << lambda << "!" << endl;
		return mu;
	}

	double v = sampling_from_gaussian(0, 1);
	double y = v * v;
	double x = mu + mu * mu * y / (2.0 *  lambda) - mu / (lambda * 2.0) * sqrt(4 * mu * lambda * y + mu * mu * y * y);

	double z = sampling_from_uniform_0_1();
	if(z <= mu / (mu + x))
		return x;
	else
		return (mu * mu) / x;
}

int sampling_from_multivector(vector<double> prob)
{
	if(prob.size() <= 0)
		return 0;
	double sum = 0;
	for(int i = 0; i < prob.size(); i++)
	{
		if(prob[i] < -1e-6)
			return 0;
		else
			sum += prob[i];
	}
	
	if(sum < 1e-6)
		return 0;
	
	for(int i = 0; i < prob.size(); i++)
		prob[i] /= sum;

	double rand_num = sampling_from_uniform_0_1();
	int cur = 0;
	while(cur < prob.size() - 1 && rand_num - prob[cur] > 0)
	{
		rand_num -= prob[cur];
		cur++;
	}
	return cur;
}
