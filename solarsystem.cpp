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

#include <stdio.h>
#include <time.h>

using std::fopen;
using std::clock;
using std::clock_t;

//Descriptive help function, explaining the basics of how to use everything.  Exits after running
void PrintHelp(void) {
    
    cout << "This is the general help section for the SolarSystem program." << endl << endl;
    cout << "Input options are as follows:" << endl;
    cout << "-h or --help: Prints this help section, overriding other options" << endl;
    cout << "filename : Only required argument.  Input file containing objects positions and " << endl;
    cout << "          velocities or orbital parameters." << endl;
    cout << "-f : Optional.  Input file format.  0=position/velocity, 1 = orbital parameters."<< endl;
    cout << "     Default is 0" << endl;
    cout << "--norecenter : Cancels recentering that is normally done." << endl;
    cout << "               Using this option may result in the system having an overall" << endl;
    cout << "               velocity drift." << endl;
    cout << "-o : Optional.  Output directory.  Default is current working directory." << endl;
    cout << "-r : Optional.  Total time to run the simulation, in years." << endl;
    cout << "     Default is 10 Earth-years." << endl;
    cout << "-s : Optional.  How many time steps between printouts.  Default is 1 (every step)" << endl;
    cout << "-t : Optional.  Size of time step (in seconds).  Default is 3600 (one hour)" << endl;
    cout << endl;
    cout << "Typical usage pattern looks like:" << endl;
    cout << "./SolarSystem filename [-f 0] [-t 3600] [-r 10] [-s 1] [-o directory]" << endl;
    
    exit(0);
}

int main(int argc, char **argv) {
    //Time to figure out the input arguments
    if(argc < 2) {
        cout << "ERROR: Insufficient number of arguments" << endl;
        cout << "FORMAT: " << endl;
        cout << "./SolarSystem filename [-f 0] [-t 3600] [-r 10] [-s 1] [-o directory]" << endl;
        cout << "If you are confused, run with the option -h to get help" << endl;
        exit(1);
    }
    
    //Initializing assorted variables
    int my_arg = 1;
    int filetype = 0;
    long double dt = 3600.;
    long double run_time = 10;
    int steps_between_prints = 1;
    string output_directory = "";
    bool do_recenter = 1;
    string filename;
    
    string temp_string;
    //Arguments needed:
    while(my_arg < argc) {
        temp_string = argv[my_arg];
        if(temp_string.find("-")==0) {
            if(temp_string.find("-f")==0) {
                //Type of position input -- X, Y, Z, or orbital elements?
                //Default X, Y, Z, filetype=0
                //Orbital elements, filetype=1
                filetype = atoi(argv[my_arg+1]);
                my_arg+=2;
            } else if(temp_string.find("-t")==0) {
                //Size of timestep (in seconds)
                //Default 1 hour
                dt = atof(argv[my_arg+1]);
                my_arg+=2;
            } else if(temp_string.find("-r")==0) {
                //Total desired runtime
                //Default is 10 Earth-years
                run_time = atof(argv[my_arg+1]);
                my_arg += 2;
            } else if(temp_string.find("-s")==0) {
                //How often to output positions
                //Default is every time step
                steps_between_prints = atof(argv[my_arg+1]);
                my_arg +=2;
            } else if(temp_string.find("-o")==0) {
                //Where to place the output files
                //Default is current working directory
                output_directory = argv[my_arg+1];
                my_arg += 2;
            } else if(temp_string.find("--norecenter")==0) {
                //Turning off recentering
                do_recenter = 0;
            } else if(temp_string.find("-h")==0 || temp_string.find("--h")==0) {
                //Catching help requests overrides any run routines -- print help and exit
                PrintHelp();
            } else {
                //Failure mode -- Doesn't understand the option.  Complain, print help, and exit
                cout << "I'm sorry, I don't understand the option " << argv[my_arg] << " that you entered.  Printing help section..." << endl;
                PrintHelp();
            }
    
        } else {
            if(filename.length() > 0) {
                cout << "WARNING: Only the first input filename will be used" << endl;
            } else {
                //Input file name
                filename = argv[1];
            }
            my_arg++;
        }
    }
    
    //Error catching for input parameter values -- revert to defaults if something's weird
    //Print a warning if we make such a change
    //Fail if there's no input file
    if(filename.length()==0) {
        cout << "ERROR!  Cannot run without an input file containing orbit information." << endl;
        exit(1);
    }
    if (dt <= 0) {
        cout << "WARNING: Timestep must be positive; reverting to default (1 hour)" << endl;
        dt = 3600.;
    }
    if (run_time <= 0) {
        cout << "WARNING: Run time must be positive.  Reverting to default (10 years)" << endl;
        run_time = 10.;
    }
    if (steps_between_prints <=0) {
        cout << "WARNING: Number of steps per printout must be at least 1.  Setting to 1." << endl;
        steps_between_prints = 1;
    }
    //And complain if we can't find the requested output directory
    FILE * testFile;
    char * output_dir_char = (char *) output_directory.c_str();
    testFile = fopen(output_dir_char, "r");
    if(testFile==NULL) {
        cout << "WARNING: Requested output directory does not exist.  Defaulting to current" << endl;
        cout << "         working directory." << endl;
        output_directory = "";
    }
    fclose(testFile);

    //Read the input file
    vector<valarray<long double> > positions;
    vector<string> names;
    cout << "Ready for readin" << endl;
    positions = ReadParameterFile(filename, filetype, names);
    cout << "Readin complete" << endl;
    
    //Initialize the solar system
    //Note that this includes setting the timestep
    SolarSystem system;
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
    if(do_recenter) system.MoveToCenterOfMass();
    
    //Run everything!  With timekeeping
    cout << "Ready to run!" << endl;
    clock_t mytime;
    mytime = clock();
    system.PrintPlanets();
    RunTheSystem(system, run_time*365.24*24.*3600.,steps_between_prints);
    mytime = mytime - clock();
    cout << "All done!  Total time used was " << mytime/CLOCKS_PER_SEC << " s." << endl;

return(0);
}
