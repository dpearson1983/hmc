#ifndef _HMC_HPP_
#define _HMC_HPP_

template <typename T> class likelihood{
    std::vector<T> x;
    std::vector<double> data, model;
    
    void calculateModel(std::vector<double> &theta);
    
    double getChisq();
    
    public:
        likelihood();
        
        void init(
        
        double getLikelihood(std::vector<double> &theta);
        
        double getNegLogLike(std::vector<double> &theta);
        
};

/* This function is to be defined by the user. Here you're model should be able to be calculated soley
 * on based on the supplied parameter vector. If you need to, you should add more data members and functions
 * in order to accomplish this calculation. For example, if you need to create a cubic spline you can add
 * new data members to store that information, and modify the constructor to initialize that object.
 */
void likelihood::calculateModel(std::vector<double> &theta) {
    
}

#endif
