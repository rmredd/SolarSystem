/*
 *  solarsystem.cpp
 *  SolarSystem
 *
 *  Created by Jadzia on 1/1/09.
 *  Copyright 2009 Rachel Reddick. All rights reserved.
 *
 */

//#include "solarsystem.h"
//#include "rungekutta.h"
//#include "predictorcorrector.h"
#include "orbital.h"

int main(int argc, char **argv) {

    //valarray<double> test = orbitalElementConverter(1.5e11,1.5e11,0,0,0,0);
    //cout<<"Orbital elements test"<<endl;
    //cout<<test[0]<<" "<<test[1]<<" "<<test[2]<<" "<<test[3]<<" "<<test[4]<<" "<<test[5]<<endl;

    //Initialize the solar system
    //Note that this includes setting the timestep
    SolarSystem system;
    double dt = 3600.;
    system.Initialize(dt, "");
    
    //Add the Sun as a planet
    valarray<double> solar_position(6);
    solar_position[0] = 0.; solar_position[1] = 0.; solar_position[2] = 0.;
    solar_position[3] = 0.; solar_position[4] = -0.0894690; solar_position[5] = 0.;
    system.AddPlanet(1.9891e30, 6.96e8, solar_position, 0., "TheSun");
    
    //Add the Earth as a planet
    valarray<double> initial(6);
    double eccEarth = 0.;//.01671;
    double semimajorEarth = 1.49498e11; //m
    double semiminorEarth = semimajorEarth*sqrt(1-eccEarth*eccEarth); //m
    double incEarth = 0; //1.+34/60. ; //degrees
    double ascNodeEarth = 0.; //348.74; //degrees
    double argPeriEarth = 0.;// 114.21; //degrees
    double startAngle = 0;
    //initial conditions set to ascending node
    initial = orbitalElementConverter(semimajorEarth,semiminorEarth,incEarth,ascNodeEarth,argPeriEarth,startAngle);
    system.AddPlanet(5.972e24, 6.371e6, initial, 0., "TheEarth");

    cout<<initial[0]<<" "<<initial[1]<<" "<<initial[2]<<" "<<initial[3]<<" "<<initial[4]<<" "<<initial[5]<<endl;
    cout<<"Radius: "<<sqrt(initial[0]*initial[0]+initial[1]*initial[1])<<endl;
    cout<<"Velocity: "<<sqrt(initial[3]*initial[3]+initial[4]*initial[4])<<endl;
    cout<<incEarth*M_PI/180<<endl;

    //Run everything!
    cout << "Ready to run!" << endl;
    system.PrintPlanets();
    RunTheSystem(system, 3*365.24*24.*3600.);
    //RunTheSystem(system, 50*3600.);

return(0);
}
