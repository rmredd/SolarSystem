//
//  predictorcorrector.h
//  SolarSystem
//
//  Created by Rachel Reddick on 2/15/14.
//
//

/*
    Predictor-corrector integration scheme, fourth order
 */

#include <valarray>

using std::valarray;

//Adams-Bashforth 4th order multistep integrator
valarray<double> AdamsBashforth(const valarray<double> xold, const double t, const valarray<double> xm1, const valarray<double> xm2, const valarray<double> xm3,double h, valarray<double>(&f)(const valarray<double>, const double)){
    const int size = xold.size();
    valarray<double> xnew(size);
    
    xnew = xold + h*( 55./24.*f(xold,t) - 59./24.*f(xm1,t-h) + 37./24.*f(xm2,t-2*h) - 0.125*f(xm3,t-3*h) );
    
    return(xnew);
}

//Adams-Moulton 4th order multistep integrator
valarray<double> AdamsMoulton(const valarray<double> xold, const double t, const valarray<double> xm1, const valarray<double> xm2, const valarray<double> xp1, double h, valarray<double>(&f)(const valarray<double>, const double)){
    const int size = xold.size();
    valarray<double> xnew(size);
    
    xnew = xold + h*( 0.125*f(xp1,t+h) + 19./24.*f(xold,t) - 5./24.*f(xm1,t-h) + 1./24.*f(xm2,t-2*h) );
    
    return(xnew);
}

//Combining the two different integrators into a single predictor-corrector
valarray<double> PredictorCorrector(const valarray<double> xold, const double t, const valarray<double> xm1, const valarray<double> xm2, const valarray<double> xm3,double h, valarray<double>(&f)(const valarray<double>, const double)){
    const int size = xold.size();
    
    valarray<double> der_xm1(size), der_xm2(size), der_xm3(size), der_xold(size), der_xp1(size), xp1, xnew;
    
    der_xold = f(xold,t);
    der_xm1 = f(xm1,t-h);
    der_xm2 = f(xm2,t-2*h);
    der_xm3 = f(xm3,t-3*h);
    
    //Adams-Bashford Prediction step
    xp1 = xold + h*(55./24.*der_xold - 59./24*der_xm1 + 37./24.*der_xm2 - 0.125*der_xm3);
    
    //Evaluation step
    der_xp1 = f(xp1,t+h);
    
    //Adams-Moulton Correction step
    xnew = xold + h*( 0.125*der_xp1 + 19./24.*der_xold - 5./24.*der_xm1 + 1./24.*der_xm2);
    
    return(xnew);
}