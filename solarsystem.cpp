/*
 *  solarsystem.cpp
 *  SolarSystem
 *
 *  Created by Jadzia on 1/1/09.
 *  Copyright 2009 Rachel Reddick. All rights reserved.
 *
 */

#include "solarsystem.h"
#include "orbital.h"
#include "parameter_readin.h"

int main(int argc, char **argv) {
    //Time to figure out the input arguments
    if(argc < 2) {
        cout << "ERROR: Insufficient number of arguments" << endl;
        cout << "FORMAT: " << endl;
        exit(1);
    }
    //Arguments needed:
    //Input file name
    string filename = argv[1];
    
    //Type of position input -- X, Y, Z, or orbital elements?
    //Default X, Y, Z
    
    //Size of timestep (in seconds)
    //Default 1 hour

    //Total desired runtime
    //Default is 3 Earth-years
    
    //How often to output positions
    //Default is every time step

    //Read the input file
    vector<valarray<long double> > positions;
    vector<string> names;
    positions = ReadParameterFile(filename, 0, names);
    
    //Initialize the solar system
    //Note that this includes setting the timestep
    SolarSystem system;
    long double dt = 3600.;
    system.Initialize(dt, "");
    
    //Add the Sun as a planet
    valarray<long double> solar_position(6);
    solar_position[0] = 0.; solar_position[1] = 0.; solar_position[2] = 0.;
    solar_position[3] = 0.; solar_position[4] = -0.0894690; solar_position[5] = 0.;
    system.AddPlanet(1.9891e30, 6.96e8, solar_position, 0., "TheSun");
    
    //Add the Earth as a planet
    valarray<long double> initial(6);
    long double eccEarth = 0.01671;
    long double semimajorEarth = 1.49498e11; //m
    long double semiminorEarth = semimajorEarth*sqrt(1-eccEarth*eccEarth); //m
    long double incEarth = 1.+34/60. ; //degrees
    long double ascNodeEarth = 348.74; //degrees
    long double argPeriEarth = 114.21; //degrees
    long double startAngle = 0;
    //initial conditions set to ascending node
    initial = orbitalElementConverter(semimajorEarth,semiminorEarth,incEarth,ascNodeEarth,argPeriEarth,startAngle);
    system.AddPlanet(5.972e24, 6.371e6, initial, 0., "TheEarth");

    cout<<initial[0]<<" "<<initial[1]<<" "<<initial[2]<<" "<<initial[3]<<" "<<initial[4]<<" "<<initial[5]<<endl;
    cout<<"Radius: "<<sqrt(initial[0]*initial[0]+initial[1]*initial[1])<<endl;
    cout<<"Velocity: "<<sqrt(initial[3]*initial[3]+initial[4]*initial[4])<<endl;
    cout<<incEarth*M_PI/180<<endl;
    
    //Add Jupiter as a planet
    initial = orbitalElementConverter(7.785e11, 7.785e11*sqrt(1-0.0488*0.0488), 1.305, 100.492, 275.066, 0);
    system.AddPlanet(1.899e27, 7.15e7, initial, 0, "Jupiter");

    //Recenter our system on the center of mass
    system.MoveToCenterOfMass();
    
    //Run everything!
    cout << "Ready to run!" << endl;
    system.PrintPlanets();
    RunTheSystem(system, 10*365.24*24.*3600.);
    //RunTheSystem(system, 50*3600.);

return(0);
}
