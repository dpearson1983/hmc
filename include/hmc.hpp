#ifndef _HMC_HPP_
#define _HMC_HPP_

#include <likelihood.hpp>

template <typename T> class hmc{
    std::vector<std::vector<double>> theta_0, theta_i, vel, acc, kick;
    std::vector<double> chisq_0, chisq_i;
    std::vector<int> trials;
    double dt;
    int steps, walkers;
    std::mt19937_64 gen;
    std::normal_distribution<double> velDist;
    std::uniform_real_distribution<double> trialDist;
    likelihood<T> L;
    
    void kicker(int tid);
    
    void updateVel(double t, int tid);
    
    void updatePos(double t, int tid);
    
    void updateACC(int tid);
    
    void integrate(int tid);
    
    bool trial();
    
    public:
        hmc();
        
        hmc(std::vector<double> theta_0, double dt, int steps, std::vector<double> &kick, likelihood<T> &L);
        
        void init(std::vector<double> theta_0, double dt, int steps, std::vector<double> & kick, 
                  likelihood<T> &L);
        
        void run(int numReals);
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
}

bool hmc<T>::trial(int tid) {
    this->trials[tid]++;
    bool moved = false;
    hmc<T>::integrate(tid);
    double L = exp(0.5*(this->chisq_0[tid] - this->chisq_i[tid]);
    double R = this->trial_dist(this->gen);
    
    if (L > R) {
        this->theta_0[tid] = this->theta_i[tid];
        this->chisq_0[tid] = this->chisq_i[tid];
    }
    return moved;
}

#endif
        
