#ifndef _HMC_HPP_
#define _HMC_HPP_

#include <likelihood.hpp>

template <typename T> class hmc{
    std::vector<double> theta_0, theta_i, vel;
    double chisq_0, chisq_i;
    double dt;
    int steps;
    likelihood<T> L;
    
    void integrate();
    
    public:
        hmc(std::vector<double> theta_0, double dt, int steps, double kick);
        
        void run(int numReals);
};

#endif
        
