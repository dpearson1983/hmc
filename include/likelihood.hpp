#ifndef _HMC_HPP_
#define _HMC_HPP_

#include <vector>
#include <string>
#include <cmath>

template <typename T> class likelihood{
    // Default data members
    std::vector<T> x;
    std::vector<double> data;
    std::vector<std::vector<double>> model;
    std::vector<std::vector<double>> Psi;
    int walkers;
    
    // User defined data members
    
    
    // Private member functions, three need to be dilled in by the user. If additional functions are needed
    // to facilitate the calculation of the model, those should be defined here, and implemented below.
    void calculateModel(std::vector<double> &theta, int tid = 0); // User defined
    
    void readDataFile(std::string file); // User defined
    
    void readCovarFile(std::string file); // User defined
    
    double getChisq(int tid = 0);
    
    public:
        likelihood();
        
        likelihood(std::string dataFile, std::string covarFile, int walkers);
        
        void init(std::string dataFile, std::string covarFile, int walkers);
        
        double getLikelihood(std::vector<double> &theta, int tid = 0);
        
        double getNegLogLike(std::vector<double> &theta, int tid = 0);
        
};

/* This function is to be defined by the user. Here you're model should be able to be calculated soley
 * on based on the supplied parameter vector. If you need to, you should add more data members and functions
 * in order to accomplish this calculation. For example, if you need to create a cubic spline you can add
 * new data members to store that information, and modify the constructor to initialize that object.
 */
void likelihood<T>::calculateModel(std::vector<double> &theta, int tid) {
    
}

/* This function should read in the data to fit the model too. It should approiately store the x values
 * and the data in the built in data members. The templated type allows you to use custom structs for
 * the x values in case you have more than one independent variable.
 */
void likelihood<T>::readDataFile(std::string file) {
    
}

/* This function should read in the covariance matrix for your data, invert it and store it into the Psi data
 * member.
 */
void likelihood<T>::readCovarFile(std::string file) {
    
}

double likelihood<T>::getChisq() {
    double chisq = 0.0;
    for (size_t i = 0; i < this->data.size(); ++i) {
        for (size_t j = i; i < this->data.size(); ++j) {
            chisqs += (this->data[i] - this->model[tid][i])*this->Psi[i][j]*(this->data[j] - this->model[tid][j]);
        }
    }
    return chisq;
}

likelihood<T>::likelihood() {
    this->walkers = 0;
}

likelihood<T>::likelihood(std::string dataFile, std::string covarFile, int walkers) {
    likelihood<T>::init(dataFile, covarFile, int walkers);
}

void likelihood<T>::init(std::string dataFile, std::string covarFile, int walkers) {
    likelihood<T>::readDataFile(dataFile);
    likelihood<T>::readCovarFile(covarFile);
    
    this->walkers = walkers;
    this->model = std::vector<std::vector<double>>(walkers, std::vector<double>(this->data.size()));
}

double likelihood<T>::getLikelihood(std::vector<double> &theta, int tid) {
    likelihood<T>::calculateModel(theta, tid);
    double chisq = likelihood<T>::getChisq(tid);
    return std::exp(-0.5*chisq);
}

double likelihood<T>::getNegLogLike(std::vector<double> &theta, int tid) {
    likelihood<T>::calculateModel(theta, tid);
    return likelihood<T>::getChisq(tid);
}

#endif
