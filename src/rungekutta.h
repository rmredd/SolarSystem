/*
 *  rungekutta.h
 *  SolarSystem
 *
 *  Created by Jadzia on 8/5/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include <valarray>
#include <vector>

using std::valarray;
using std::vector;

//Runge Kutta integration step
valarray<double> RungeKutta(const valarray<double> xold, const double t, const double h, valarray<double>(&f)(const valarray<double>, const double)) {
  const int size = xold.size();
  valarray<double> k1(size), k2(size), k3(size), k4(size);
  k1=h*f(xold,t);
  k2=h*f(xold+k1,t+h/2);
  k3=h*f(xold+(.5)*k2,t+h/2);
  k4=h*f(xold+(.5)*k3,t+h);
  valarray<double> xnew(size);
  xnew = xold+(k1+(2.0)*k2+(2.0)*k3+k4)/6;
  return(xnew);
}