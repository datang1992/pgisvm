

#ifndef __STORED__INFO___
#define __STORED__INFO___



#include <vector>
#include <math.h>
#include "matrix.h"
#include <gsl/gsl_randist.h>

using namespace std;



class stored_info
{
	public:
		int K, N;

		int dim;

		vector<double> y;

		vector<matrix> x;

		vector<int> z, z_l, z_r, z_s;
		vector<double> zeta;
		vector<int> Num, Num_l, Num_r;
		vector<double> pi, pi_l, pi_r;		
		vector<int> lr;		
		//vector<double> split_probability;
		vector<matrix> eta, eta_l, eta_r;
		vector<matrix> gamma_mu, gamma_sigma, gamma_mu_l, gamma_sigma_l, gamma_mu_r, gamma_sigma_r;

		vector<matrix> gamma_mu0, gamma_mu0_l, gamma_mu0_r, gamma_mu0_s;
		vector<matrix> gamma_psi, gamma_psi_l, gamma_psi_r, gamma_psi_s;
		vector<double> gamma_m, gamma_m_l, gamma_m_r, gamma_m_s;
		vector<double> gamma_n0, gamma_n0_l, gamma_n0_r, gamma_n0_s;

		int initial_clusters_number;
		

		matrix initial_mu0, initial_psi;
		double initial_m, initial_n0;
		
		matrix mark;

		vector<double> split_probability;

		int Nthreads;
		vector<double> splittable;
		
		int number_of_iterations;


		bool useSuperclusters;
		vector<int> superclusters, supercluster_labels;

		vector<double> omega;

		double C, l, nu;
		double alpha;
		
		long rand_seed;
		gsl_rng	*rng;

		

		void initialize();
		void sample_cluster();
		void sample_superclusters();
		void sample_label();
		void propose_merges();
		void propose_splits();
		void propose_random_split_assignments();
		

        stored_info();
		stored_info(int, int);
		~stored_info();
};


		













#endif
