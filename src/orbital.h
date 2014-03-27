/*
 *  orbital.h
 *  SolarSystem
 *
 *  Created by Jadzia on 8/5/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include <valarray>
#include <iostream>

using std::valarray;
using std::cout;
using std::endl;

//Making sure that G and Msun are defined
#ifndef G
#define G 6.67384e-11 //m^3/(kg s^2)
#endif

#ifndef Msun
#define Msun 1.989e30 //kg
#endif

//Function prototypes
valarray<long double> gravitySun(const valarray<long double>, const long double);

void rotateZ(long double& x, long double& y, long double& z, long double theta);
void rotateX(long double& x, long double& y, long double& z, long double theta);
valarray<long double> orbitalElementConverter(long double, long double, long double, long double, long double, long double);