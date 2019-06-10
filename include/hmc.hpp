#ifndef _HMC_HPP_
#define _HMC_HPP_

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <random>
#include <cmath>
#include <likelihood.hpp>

template <typename T> class hmc{
    std::vector<std::vector<double>> theta_0, theta_i, vel, acc, kick;
    std::vector<double> chisq_0, chisq_i;
    std::vector<int> trials, accepted;
    double dt;
    int steps, walkers, digits;
    std::mt19937_64 gen;
    std::normal_distribution<double> velDist;
    std::uniform_real_distribution<double> trialDist;
    likelihood<T> L;
    
    void kicker(int tid);
    
    void updateVel(double t, int tid);
    
    void updatePos(double t, int tid);
    
    void updateACC(int tid);
    
    void integrate(int tid);
    
    void writeThetaScreen(int tid)
    
    void trial();
    
    int numDigits();
    
    std::string filename(std::string base, int tid, std::string ext);
    
    public:
        hmc();
        
        hmc(std::vector<double> &theta_0, double dt, int steps, std::vector<double> &kick, likelihood<T> &L);
        
        void init(std::vector<double> &theta_0, double dt, int steps, std::vector<double> & kick, 
                  likelihood<T> &L);
        
        void run(int numReals, int numBurn);
};

void hmc<T>::kicker(int tid) {
    for (size_t i = 0; i < this->theta_0[tid].size(); ++i) {
        this->vel[tid][i] = this->kick[i]*this->velDist(this->gen);
    }
}

void hmc<T>::updateVel(double t, int tid) {
    for (size_t i = 0; i < this->theta_i[tid].size(); ++i) {
        this->vel[tid][i] += 0.5*this->acc[tid][i]*t;
    }
}

void hmc<T>::updatePos(double t, int tid) {
    for (size_t i = 0; i < this->theta_i[tid].size(); ++i) {
        this->theta_i[tid][i] += this->vel[tid][i]*t + 0.5*this->acc[tid][i]*t;
    }
}

void hmc<T>::updateAcc(int tid) {
    double h = 1E-6;
    for (size_t i = 0; i < this->theta_i[tid].size(); ++i) {
        this->theta_i[tid][i] += 0.5*h;
        double f_p = this->L.getNegLogLike(this->theta_i[tid], tid);
        this->theta_i[tid][i] -= h;
        double f_m = this->L.getNegLogLike(this->theta_i[tid], tid);
        this->theta_i[tid][i] += 0.5*h;
        this->acc[tid][i] = (f_p - f_m)/h;
    }
}

void hmc<T>::integrate(int tid) {
    this->theta_i = this->theta_0;
    hmc<T>::kicker(tid);
    
    for (int i = 0; i < steps; ++i) {
        hmc<T>::updateVel(0.5*this->dt, tid);
        hmc<T>::updatePos(this->dt, tid);
        hmc<T>::updateAcc(tid);
        hmc<T>::updateVel(0.5*this->dt, tid);
    }
    this->chisq_i = L.getNegLogLike(this->theta_i[tid], tid);
}

void hmc<T>::writeThetaScreen(int tid) {
    if (tid == 0) {
        std::cout << "\r";
        std::cout << std::setw(15) << this->trials[tid];
        for (size_t i = 0; i < this->theta_0[tid].size(); ++i) {
            std::cout << std::setw(15) << this->theta_0[tid][i];
        }
        std::cout << std::setw(15) << this->chisq_0[tid];
    }
}

void hmc<T>::trial(int tid) {
    this->trials[tid]++;
    hmc<T>::integrate(tid);
    double L = exp(0.5*(this->chisq_0[tid] - this->chisq_i[tid]);
    double R = this->trial_dist(this->gen);
    
    if (L > R) {
        this->theta_0[tid] = this->theta_i[tid];
        this->chisq_0[tid] = this->chisq_i[tid];
        hmc<T>::writeThetaScreen(tid);
        this->accepted[tid]++;
    }
}

void hmc<T>::numDigits() {
    this->digits = 0;
    int temp = this->walkers;
    while(temp) {
        temp /= 10;
        this->digits++;
    }
}

std::string hmc<T>::filename(std::string base, int tid, std::string ext) {
    std::stringstream file;
    file << base << std::setw(this->digits) << std::setfill('0') << tid << "." << ext;
    return file.str();
}

// TODO: Need to call the constructors for the random number generator and the distributions
hmc<T>::hmc() {
    this->walkers = 0;
}

// TODO: Need to call the constructors for the random number generator and the distributions
hmc<T>::hmc(std::vector<double> &theta_0, double dt, int steps, std::vector<double> &kick, likelihood<T> &L) {
    hmc<T>::init(theta_0, dt, steps, kick, L);
}

hmc<T>::init(std::vector<double> &theta_0, double dt, int steps, std::vector<double> &kick, likelihood<T> &L) {
    this->walkers = walkers;
    hmc<T>::numDigits();
    this->steps = steps;
    this->dt = dt;
    this->theta_i = std::vector<std::vector<double>>(walkers, std::vector<double>(theta_0.size()));
    this->vel = std::vector<std::vector<double>>(walkers, std::vector<double>(theta_0.size()));
    this->acc = std::vector<std::vector<double>>(walkers, std::vector<double>(theta_0.size()));
    this->kick = kick;
    this->chisq_0 = std::vector<double>(walkers);
    this->chisq_i = std::vector<double>(walkers);
    this->trials = std::vector<int>(walkers);
    this->L = L;
    
    for (int i = 0; i < walkers; ++i) {
        this->theta_0.push_back(theta_0);
    }
}

void hmc<T>::run(int numReals, int numBurn, std::string chainBase, std::string chainExt) {
    omp_set_num_threads(this->walkers);
    
    std::cout << "Buring the first " << numBurn << " trials from each chain..." << std::endl;
#pragma omp parallel 
    {
        int tid = omp_get_thread_num();
        for (int i = 0; i < numBurn; ++i) {
            hmc<T>::trial(tid);
        }
    }
    
    std::cout << "Burn in acceptance rates from each walker:\n";
    for (int i = 0; i < this->walkers; ++i) {
        double acceptance = double(this->accepted[i])/double(this->trials[i]);
        std::cout << "Walker #" << i + 1 << " = " << acceptance << std::endl;
    }
    
    std::cout << "Running the chains..." << std::endl
#pragma omp parallel
    {
        this->trials[tid] = 0;
        this->accepted[tid] = 0;
        std::string chainFile = hmc<T>::filename(chainBase, tid, chainExt);
        std::ofstream fout(chainFile);
        for (int i = 0; i < numReals; ++i) {
            hmc<T>::trial(tid);
            fout << i + 1 << " ";
            for (size_t par = 0; par < this->theta_0[tid].size(); ++par)
                fout << this->theta_0[tid][par] << " ";
            fout << this->chisq_0[tid] << "\n";
        }
        fout.close();
    }
}

#endif
        
