
#include "train.h"
#include "matrix.h"
#include "stored_info.h"
#include "utilities.h"
#include "sampling.h"
#include <gsl/gsl_sf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_cdf.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <vector>
#include <iostream>
using namespace std;


void initialization(stored_info *info) {
    cout << "Initialization start!\n";
    info -> K = info -> initial_clusters_number;
    //cout << gsl_rng_uniform(info -> rng);
    matrix tmp = matrix(info -> dim, 1), identity = matrix(info -> dim, -2);
    //cout << tmp.m << ' ' << tmp.n << ' ' << identity.m << identity.n << endl;
    //cout << info -> nu;
    for (int i = 0; i < info -> K; i++) {
        info -> Num.push_back(0);
        info -> Num_l.push_back(0);
        info -> Num_r.push_back(0);
        //info -> pi.push_back(0);
        //info -> pi_l.push_back(0);
        //info -> pi_r.push_back(0);
        info -> eta.push_back(sampling_from_gaussian(tmp, identity * ((info -> nu) * (info -> nu))));
        info -> b.push_back(sampling_from_gaussian(0, (info -> nu2) * (info -> nu2)));
	info -> eta_l.push_back(tmp);
        info -> eta_r.push_back(tmp);
	info -> b_l.push_back(0);
	info -> b_r.push_back(0);
        info -> gamma_m.push_back(info -> initial_m);
        info -> gamma_mu0.push_back(info -> initial_mu0);
        info -> gamma_n0.push_back(info -> initial_n0);
        info -> gamma_psi.push_back(info -> initial_psi);
        info -> split_probability.push_back(0);
    }
    cout << "Finish pushing back variables!\n";
    
    double *alpha = new double[info -> K], *theta = new double[info -> K];
    for (int i = 0; i < info -> K; i++) {
        alpha[i] = info -> alpha;
    }
    gsl_ran_dirichlet(info -> rng, info -> K, alpha, theta);
    for (int i = 0; i < info -> K; i++) {
        info -> pi.push_back(theta[i]);
        double alpha_sub[2] = {info -> alpha, info -> alpha}, theta_sub[2];
        gsl_ran_dirichlet(info -> rng, 2, alpha_sub, theta_sub);
        info -> pi_l.push_back(theta_sub[0]);
        info -> pi_r.push_back(theta_sub[1]);
        ///cout << info -> pi[i] << " ";
    }
    
    vector<double> prob;
    
    delete[] alpha;
    delete[] theta;
    
    cout << "Finish sampling weights!\n";
    
    //cout << prob.size() << ' ' << info -> K << endl;
    
    for (int i = 0; i < info -> K; i++) {
        prob.push_back(theta[i]);
    }
    
    //#pragma omp parallel for
    for (int i = 0; i < info -> N; i++) {
        info -> z.push_back(sampling_from_multivector(prob));
        info -> omega.push_back(1);
        //cout << info -> z[i] << ' ';
        info -> Num[info -> z[i]]++;
        vector<double> prob_sub;
        prob_sub.push_back(info -> pi_l[info -> z[i]]);
        prob_sub.push_back(info -> pi_r[info -> z[i]]);
        info -> lr.push_back(1 - sampling_from_multivector(prob_sub));
        info -> zeta.push_back(1);
        if (info -> lr[i]) {
            info -> Num_l[info -> z[i]]++;
        }
        else
            info -> Num_r[info -> z[i]]++;
    }
    cout << "Initialization finished!\n";
}

void sampling_cluster_parameters(stored_info *info) {
    double PI = 3.14159265358979323846264338327950288;
    //#pragma omp parallel for
    for (int i = 0; i < info -> K; i++) {
        double split_prob = log(info -> alpha);
        matrix sigma = matrix(info -> dim, -2) * (1.0 / (info -> nu * info -> nu));
        double sigma2 = 1 / ((info -> nu2) * (info -> nu2));
        for (int j = 0; j < info -> N; j++) {
            if (info -> z[j] == i) {
                sigma = sigma + ((info -> x[j]) * (info -> x[j]).transpose()) * ((info -> C) * (info -> C) / info -> omega[j]);
                sigma2 = sigma2 + ((info -> C) * (info -> C) / info -> omega[j]);
            }
        }
        sigma = sigma.inverse();
        sigma2 = 1 / sigma2;
        matrix mu = matrix(info -> dim, 1);
        double mu2 = 0;
        for (int j = 0; j < info -> N; j++) {
            if (info -> z[j] == i) {
                mu = mu + sigma * (info -> x[j]) * ((info -> C) * (info -> y[j]) * ((info -> omega[j]) + (info -> C) * ((info -> l) - (info -> y[j]) * (info -> b[i]))) / (info -> omega[j]));
                mu2 = mu2 + sigma2 * ((info -> C) * (info -> y[j]) * ((info -> omega[j]) + (info -> C) * ((info -> l) - (info -> y[j]) * ((info -> eta[i]).transpose() * (info -> x[j])).num[0][0])) / (info -> omega[j]));
            }
        }
        info -> eta[i] = sampling_from_gaussian(mu, sigma);
        info -> b[i] = sampling_from_gaussian(mu2, sigma2);
        if (i == 0) {
            cout << mu2 << ' ' << sigma2 << ' ' << info -> b[i] << ' ' << mu.num[0][0] << ' ' << sigma.num[0][0] << endl;
        }
        
        info -> Num[i] = 0;
        info -> Num_l[i] = 0;
        info -> Num_r[i] = 0;
        for (int j = 0; j < info -> N; j++) {
            if (info -> z[j] == i) {
                info -> Num[i]++;
                if (info -> lr[j] == 1) {
                    info -> Num_l[i]++;
                }
                else
                    info -> Num_r[i]++;
            }
        }

        if (info -> Num_l[i] && info -> Num_r[i]) {
            split_prob += gsl_sf_lngamma(info -> Num_l[i]) + gsl_sf_lngamma(info -> Num_r[i]) - gsl_sf_lngamma(info -> Num[i]);
        }
        split_prob -= 0.5 * log(sigma2) + 0.5 * log(sigma.det()) - log(info -> nu2) - (info -> dim) * log(info -> nu) + 0.5 * (mu.transpose() * sigma.inverse() * mu).num[0][0] + 0.5 * (mu2 * mu2 / sigma2);
        for (int j = 0; j < info -> N; j++) {
            if (info -> z[j] == i) {
                split_prob -= -log(2 * PI * (info -> omega[j])) - 0.5 * ((info -> omega[j]) + (info -> C) * (info -> l)) * ((info -> omega[j]) + (info -> C) * (info -> l)) / (info -> omega[j]);
            }
        }
        
        double alpha[2] = {info -> Num_l[i] + info -> alpha, info -> Num_r[i] + info -> alpha};
        double theta[2];
        gsl_ran_dirichlet(info -> rng, 2, alpha, theta);
        info -> pi_l[i] = theta[0];
        info -> pi_r[i] = theta[1];
        
        matrix sigma_l = matrix(info -> dim, -2) * (1.0 / ((info -> nu) * (info -> nu)));
        double sigma2_l = 1.0 / ((info -> nu2) * (info -> nu2));
        for (int j = 0; j < info -> N; j++) {
            if (info -> z[j] == i && info -> lr[j] == 1) {
                sigma_l = sigma_l + ((info -> x[j]) * (info -> x[j]).transpose()) * ((info -> C) * (info -> C) / info -> omega[j]);
                sigma2_l = sigma2_l + ((info -> C) * (info -> C) / info -> omega[j]);
            }
        }
        sigma_l = sigma_l.inverse();
        sigma2_l = 1 / sigma2_l;
        matrix mu_l = matrix(info -> dim, 1);
        double mu2_l = 0;
        for (int j = 0; j < info -> N; j++) {
            if (info -> z[j] == i && info -> lr[j] == 1) {
                mu_l = mu_l + sigma_l * (info -> x[j]) * ((info -> C) * (info -> y[j]) * ((info -> omega[j]) + (info -> C) * ((info -> l) - (info -> y[j]) * (info -> b[i]))) / (info -> omega[j]));
                mu2_l = mu2_l + sigma2_l * ((info -> C) * (info -> y[j]) * ((info -> omega[j]) + (info -> C) * ((info -> l) - (info -> y[j]) * ((info -> eta[i]).transpose() * (info -> x[j])).num[0][0])) / (info -> omega[j]));
            }
        }
        info -> eta_l[i] = sampling_from_gaussian(mu_l, sigma_l);
        info -> b_l[i] = sampling_from_gaussian(mu2_l, sigma2_l);
        
        split_prob += 0.5 * log(sigma2_l) + 0.5 * log(sigma_l.det()) - log(info -> nu2) - (info -> dim) * log(info -> nu) + 0.5 * (mu_l.transpose() * sigma_l.inverse() * mu_l).num[0][0] + 0.5 * (mu2_l * mu2_l / sigma2_l);
        for (int j = 0; j < info -> N; j++) {
            if (info -> z[j] == i && info -> lr[j] == 1) {
                split_prob += -log(2 * PI * (info -> omega[j])) - 0.5 * ((info -> omega[j]) + (info -> C) * (info -> l)) * ((info -> omega[j]) + (info -> C) * (info -> l)) / (info -> omega[j]);
            }
        }
        
        matrix sigma_r = matrix(info -> dim, -2) * (1.0 / (info -> nu * info -> nu));
        double sigma2_r = 1 / ((info -> nu2) * (info -> nu2));
        for (int j = 0; j < info -> N; j++) {
            if (info -> z[j] == i && info -> lr[j] == 0) {
                sigma_r = sigma_r + ((info -> x[j]) * (info -> x[j]).transpose()) * ((info -> C) * (info -> C) / info -> omega[j]);
                sigma2_r = sigma2_r + ((info -> C) * (info -> C) / info -> omega[j]);
            }
        }
        sigma_r = sigma_r.inverse();
        sigma2_r = 1 / sigma2_r;
        matrix mu_r = matrix(info -> dim, 1);
        double mu2_r = 0;
        for (int j = 0; j < info -> N; j++) {
            if (info -> z[j] == i && info -> lr[j] == 0) {
                mu_r = mu_r + sigma_r * (info -> x[j]) * ((info -> C) * (info -> y[j]) * ((info -> omega[j]) + (info -> C) * ((info -> l) - (info -> y[j]) * (info -> b[i]))) / (info -> omega[j]));
                mu2_r = mu2_r + sigma2_r * ((info -> C) * (info -> y[j]) * ((info -> omega[j]) + (info -> C) * ((info -> l) - (info -> y[j]) * ((info -> eta[i]).transpose() * (info -> x[j])).num[0][0])) / (info -> omega[j]));
            }
        }
        info -> eta_r[i] = sampling_from_gaussian(mu_r, sigma_r);
        info -> b_r[i] = sampling_from_gaussian(mu2_r, sigma2_r);
        
        split_prob += 0.5 * log(sigma2_r) + 0.5 * log(sigma_r.det()) - log(info -> nu2) - (info -> dim) * log(info -> nu) + 0.5 * (mu_r.transpose() * sigma_r.inverse() * mu_r).num[0][0] + 0.5 * (mu2_r * mu2_r / sigma2_r);
        for (int j = 0; j < info -> N; j++) {
            if (info -> z[j] == i && info -> lr[j] == 0) {
                split_prob += -log(2 * PI * (info -> omega[j])) - 0.5 * ((info -> omega[j]) + (info -> C) * (info -> l)) * ((info -> omega[j]) + (info -> C) * (info -> l)) / (info -> omega[j]);
            }
        }
        
        if (exp(split_prob) < 1) {
            info -> split_probability[i] = exp(split_prob);
        }
        else
            info -> split_probability[i] = 1;
        /*
        vector<matrix> X;
        for (int j = 0; j < info -> N; j++) {
            if (info -> z[j] == i) {
                X.push_back(info -> x[j]);
            }
        }
        
        cal_niw_posterior(info -> initial_mu0, info -> initial_m, info -> initial_psi, info -> initial_n0, X, info -> gamma_mu0[i], info -> gamma_m[i], info -> gamma_psi[i], info -> gamma_n0[i]);
        */
        
    }
}

void DFS(int now, matrix &new_mark, int K, bool *Mark, int *node, int &p) {
    bool finish = true;
    int next;
    for (int i = 0; i < K; i++) {
        if (!Mark[i] && new_mark.num[now][i]) {
            finish = false;
            next = i;
            break;
        }
    }
    if (finish) {
        for (int i = 0; i < p; i++) {
            for (int j = 0; j < p; j++) {
                new_mark.num[i][j] = 1;
            }
        }
    }
    else {
        Mark[next] = true;
        node[p++] = next;
        DFS(next, new_mark, K, Mark, node, p);
    }
}

void sampling_supercluster_parameters(stored_info *info) {
    matrix distance = matrix(info -> K, info -> K);
    #pragma omp parallel for
    for (int i = 0; i < (info -> K) * (info -> K); i++) {
        int i1 = i % (info -> K), i2 = i / (info -> K);
        distance.num[i1][i2] = 1 / (((info -> gamma_mu0[i1]) - (info -> gamma_mu0[i2])).transpose() * ((info -> gamma_mu0[i1]) - (info -> gamma_mu0[i2]))).num[0][0];
    }
    matrix new_mark = matrix(info -> K, info -> K);
    #pragma omp parallel for
    for (int i = 0; i < info -> K; i++) {
        vector<double> prob;
        for (int j = 0; j < info -> K; j++) {
            prob.push_back(distance.num[i][j]);
        }
        int neighbor = sampling_from_multivector(prob);
        new_mark.num[i][neighbor] = new_mark.num[neighbor][i] = 1;
    }
    
    bool *Mark = new bool[info -> K];
    int *node = new int[info -> K], p;
    #pragma omp parallel for
    for (int i = 0; i < info -> K; i++) {
        Mark[i] = false;
    }
    for (int i = 0; i < info -> K; i++) {
        if (!Mark[i]) {
            Mark[i] = true;
            p = 0;
            DFS(i, new_mark, info -> K, Mark, node, p);
        }
    }
    delete[] Mark;
    delete[] node;
    
    info -> mark = new_mark;
}

void sampling_data_parameters(stored_info *info) {
    double PI = 3.14159265358979323846264338327950288;
    #pragma omp parallel for
    for (int i = 0; i < info -> N; i++) {
        vector<double> prob;
        vector<int> index;
        for (int j = 0; j < info -> K; j++) {
            if (info -> mark.num[info -> z[i]][j]) {
                double zeta = (info -> l) - (info -> y[i]) * (((info -> eta[j]).transpose() * (info -> x[i])).num[0][0] + (info -> b[j]));
                prob.push_back((info -> pi[j]) / (sqrt(2 * PI * info -> omega[i])) * exp(-(info -> omega[i] + (info -> C) * zeta) * (info -> omega[i] + (info -> C) * zeta) / (2 * info -> omega[i])));
                index.push_back(j);
            }
        }
        info -> z[i] = index[sampling_from_multivector(prob)];
        //cout << info -> eta.size() << ' ' << info -> z[i] << endl;
        
        info -> zeta[i] = (info -> l) - (info -> y[i]) * (((info -> eta[info -> z[i]]).transpose() * (info -> x[i])).num[0][0] + (info -> b[info -> z[i]]));
        prob.clear();
        
        double zeta = (info -> l) - (info -> y[i]) * (((info -> eta_l[info -> z[i]]).transpose() * (info -> x[i])).num[0][0] + (info -> b_l[info -> z[i]]));
        prob.push_back((info -> pi_l[i]) / (sqrt(2 * PI * info -> omega[i])) * exp(-(info -> omega[i] + (info -> C) * zeta) * (info -> omega[i] + (info -> C) * zeta) / (2 * info -> omega[i])));
        prob.clear();
        zeta = (info -> l) - (info -> y[i]) * (((info -> eta_r[info -> z[i]]).transpose() * (info -> x[i])).num[0][0] + (info -> b_r[info -> z[i]]));
        prob.push_back((info -> pi_r[i]) / (sqrt(2 * PI * info -> omega[i])) * exp(-(info -> omega[i] + (info -> C) * zeta) * (info -> omega[i] + (info -> C) * zeta) / (2 * info -> omega[i])));
        info -> lr[i] = 1 - sampling_from_multivector(prob);
    }
    
    #pragma omp parallel for
    for (int i = 0; i < info -> N; i++) {
        info -> omega[i] = 1 / sampling_from_inverse_gaussian(1 / (info -> C) / fabs(info -> zeta[i]), 1);
    }
    
}

void splitting_clusters(stored_info *info) {
    #pragma omp parallel for
    for (int i = 0; i < info -> K; i++) {
        if (gsl_rng_uniform(info -> rng) < (info -> split_probability[i])) {
            int i2;
            #pragma omp critical
            {
                matrix tmp = matrix(info -> dim, 1);
                i2 = info -> K;
                //cout << info -> K << endl;
                info -> Num.push_back(info -> Num_r[i]);
                info -> Num_l.push_back(0);
                info -> Num_r.push_back(0);
                info -> pi.push_back(0);
                info -> pi_l.push_back(0);
                info -> pi_r.push_back(0);
                info -> eta.push_back(info -> eta_r[i]);
                info -> eta_l.push_back(tmp);
                info -> eta_r.push_back(tmp);
                info -> b.push_back(info -> b_r[i]);
                info -> b_l.push_back(0);
                info -> b_r.push_back(0);
                info -> gamma_m.push_back(info -> initial_m);
                info -> gamma_mu0.push_back(info -> initial_mu0);
                info -> gamma_n0.push_back(info -> initial_n0);
                info -> gamma_psi.push_back(info -> initial_psi);
                info -> split_probability.push_back(0);
            }
            info -> Num[i] = info -> Num_l[i];
            info -> eta[i] = info -> eta_l[i];
            double alpha[2] = {info -> alpha, info -> alpha};
            double theta_l[2], theta_r[2];
            gsl_ran_dirichlet(info -> rng, 2, alpha, theta_l);
            gsl_ran_dirichlet(info -> rng, 2, alpha, theta_r);
            info -> pi_l[i] = theta_l[0];
            info -> pi_r[i] = theta_l[1];
            info -> pi_l[i2] = theta_r[0];
            info -> pi_r[i2] = theta_r[1];
            for (int j = 0; j < info -> N; j++) {
                if (info -> z[j] == i) {
                    vector<double> prob;
                    if (info -> lr[j] == 1) {
                        prob.push_back(info -> pi_l[i]);
                        prob.push_back(info -> pi_r[i]);
                        info -> lr[j] = 1 - sampling_from_multivector(prob);
                    }
                    else if(info -> lr[j] == 0) {
                        info -> z[j] = i2;
                        prob.push_back(info -> pi_l[i2]);
                        prob.push_back(info -> pi_r[i2]);
                        info -> lr[j] = 1 - sampling_from_multivector(prob);
                    }
                }
            }
        }
    }
}

void train(stored_info *info) {
    initialization(info);
    
    for (int i = 0; i < info -> number_of_iterations; i++) {
        //cout << "Iteration " << i << " starts!\n";
        double *alpha = new double[(info -> K) + 1], *theta = new double[(info -> K) + 1];
        for (int i = 0; i < info -> K; i++) {
            alpha[i] = info -> Num[i];
        }
        alpha[info -> K] = info -> alpha;
        gsl_ran_dirichlet(info -> rng, (info -> K) + 1, alpha, theta);
        for (int j = 0; j < info -> K; j++) {
            info -> pi[j] = theta[j];
        }
        delete[] alpha;
        delete[] theta;
        sampling_cluster_parameters(info);
        //cout << "Finish sampling cluster parameters!\n";
        sampling_supercluster_parameters(info);
        //cout << "Finish sampling supercluster parameters!\n";
        sampling_data_parameters(info);
        //cout << "Finish sampling data parameters!\n";
        splitting_clusters(info);
        cout << "Iteration " << i << " finishes!" << endl;
    }
    
}

double cross_validiction(stored_info *info) {
    train(info);
    
    double correct = 0;
    for (int i = 0; i < info -> N; i++) {
        int prediction = ((((info -> eta[info -> z[i]]).transpose() * (info -> x[i])).num[0][0] + 1 * (info -> b[info -> z[i]])) >= 0) ? 1 : 0;
        if (prediction == info -> y[i]) {
            correct += 1;
        }
    }
    return correct / info -> N;
    
    //return 0;
}


