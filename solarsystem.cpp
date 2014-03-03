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
    
    int myarg = 0;
    int filetype = 0;
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
    filetype = 0;
    
    //Where to place the output files
    //Default is current working directory
    string output_directory = "";

    //Read the input file
    vector<valarray<long double> > positions;
    vector<string> names;
    cout << "Ready for readin" << endl;
    positions = ReadParameterFile(filename, filetype, names);
    cout << "Readin complete" << endl;
    
    //Initialize the solar system
    //Note that this includes setting the timestep
    SolarSystem system;
    long double dt = 3600.;
    system.Initialize(dt, output_directory);

    //Add all the desired input objects
    valarray<long double> xtemp(6);
    for(int i=0; i<positions.size(); i++) {
        for(int j=0; j<6; j++){
            xtemp[j] = positions[i][j+2];
        }
        system.AddPlanet(positions[i][0], positions[i][1], xtemp, 0, names[i]);
    }
    cout << "Done initializing " << positions.size() << " objects" << endl;

    //Recenter our system on the center of mass
    system.MoveToCenterOfMass();
    
    //Run everything!
    cout << "Ready to run!" << endl;
    system.PrintPlanets();
    RunTheSystem(system, 10*365.24*24.*3600.);
    //RunTheSystem(system, 50*3600.);

return(0);
}
