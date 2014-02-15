/*
 *  solarsystem.cpp
 *  SolarSystem
 *
 *  Created by Jadzia on 1/1/09.
 *  Copyright 2009 Rachel Reddick. All rights reserved.
 *
 */

//#include "solarsystem.h"
#include "rungekutta.h"
#include "orbital.h"

int main(int argc, char **argv) {

valarray<double> test = orbitalElementConverter(1.5e11,1.5e11,0,0,0,0);
cout<<"Orbital elements test"<<endl;
cout<<test[0]<<" "<<test[1]<<" "<<test[2]<<" "<<test[3]<<" "<<test[4]<<" "<<test[5]<<endl;


//Initialize Earth
Planet earth;
valarray<double> initial(6);
double eccEarth = .01671;
double semimajorEarth = 1.49498e11; //m
double semiminorEarth = semimajorEarth*sqrt(1-eccEarth*eccEarth); //m
double incEarth = 1.+34/60. ; //degrees
double ascNodeEarth = 348.74; //degrees
double argPeriEarth = 114.21; //degrees
double startAngle = 0;
//initial conditions set to ascending node
initial = orbitalElementConverter(semimajorEarth,semiminorEarth,incEarth,ascNodeEarth,argPeriEarth,startAngle);
earth.Initialize(5.974e24,initial,0);

cout<<initial[0]<<" "<<initial[1]<<" "<<initial[2]<<" "<<initial[3]<<" "<<initial[4]<<" "<<initial[5]<<endl;
cout<<"Radius: "<<sqrt(initial[0]*initial[0]+initial[1]*initial[1])<<endl;
cout<<"Velocity: "<<sqrt(initial[3]*initial[3]+initial[4]*initial[4])<<endl;
cout<<incEarth*M_PI/180<<endl;

//initial test
double dt = 3600;
valarray<double> xold = initial;
valarray<double> xnew(6);
valarray<double> accel(6);
for (int i=0; i<10*365*24*3600/dt; i++) {
//for(int i=0; i<3*3600/dt; i++) { //short test loop
	//cout<<"okay "<<i<<endl;
	xnew = RungeKutta(xold,i*dt,dt,gravitySun);
    accel = gravitySun(xold, 0.);
    cout<<"TEST: "<<i<<" "<<accel[0]<<" "<<accel[1]<<" "<<accel[2]<<" "<<accel[3]<<" "<<endl;
	//cout<<xnew[5]-xold[5]<<endl;
	earth.Iterate(xnew,(i+1)*dt);
	
	xold = xnew;
}
cout<<"done iterating"<<endl;

earth.PrintOutput("test_output.dat");

return(0);
}
