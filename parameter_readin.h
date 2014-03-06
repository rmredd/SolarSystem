//
//  parameter_readin.h
//  SolarSystem
//
//  Created by Rachel Reddick on 2/22/14.
//
//

//General function that manages all our basic readin stuffs

#include <valarray>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>

#include "orbital.h"

using std::vector;
using std::valarray;
using std::fstream;
using std::ifstream;
using std::fscanf;
using std::string;

//This is the basic function that reads in a file; it reads in the first entry in each line as a name,
//and the remaining nentries as long double inputs for later
//Comments are indicated with "#"
vector<valarray<long double> > ParameterFileBasicRead(string filename, int nentries, vector<string> & names){

    //Figure out how big an object we'll need
    ifstream input;
    input.open(filename);
    string myline;
    int nlines = 0;
    while(!input.eof()) {
        getline(input,myline);
        if(myline.find("#")==0) continue;
        nlines++;
    }
    input.close();
    
    vector<valarray<long double> > positions(nlines);
    valarray<long double> xtemp(nentries);
    vector<string> mynames(nlines);
    
    int linecount = 0, place, place_end;
    input.open(filename);
    for(int i=0; i<nlines; i++) {
        getline(input,myline);
        //Extract the name of the planet first
        mynames[i] = myline.substr(0,myline.find(" "));
        myline.erase(0,myline.find(" "));
        while(myline.find(" ")==0) myline.erase(0,1);
        place=0;
        for(int j=0; j<8; j++){
            place_end = myline.find(" ",place);
            if(place_end==string::npos) place_end = myline.length()-1;
            if(j<7 && place_end == myline.length()-1) {
                cout << "WARNING: Error in file readin at line " << i << endl;
                cout << "         Insufficient entries in line for chosen format" << endl;
                break;
            }
            xtemp[j] = atof(myline.substr(place,place_end).c_str());
            place = place_end;
            while(myline.find(" ",place)==place){
                place++;
            }
        }
        positions[i] = xtemp;
    }
    input.close();
    
    //Fill out the names vector
    names = mynames;
    
    return(positions);
}

vector<valarray<long double> > ReadParameterFile(string filename, int filetype, vector<string> & names) {
    
    //Note that the positions valarray contains mass, radius, followed by position information
    vector<valarray<long double> > positions;
    switch (filetype) {
        case 0:
        {
            //Name(string), mass, radius, X, Y, Z, VX, VY, VZ position data
            positions = ParameterFileBasicRead(filename, 8, names);
        }
            break;
        case 1:
            //Orbital elements -- will need to be converted at some point
            //Input lines must contain name, mass, radius, semimajor axis, eccentricity, inclination, ascending node, argument of perihelion, angle
        {
            positions = ParameterFileBasicRead(filename, 8, names);
            //Run the conversion from orbital elements to X, Y, Z, VX, VY, VZ
            valarray<long double> xtemp(6);
            for(int i=0; i<positions.size(); i++) {
                for(int j=0; j<6; j++) xtemp[j] = positions[i][2+j];
                //Convert eccentricity to semiminor axis length
                xtemp[1] = xtemp[0]*(1-xtemp[1]*xtemp[1]);
                //Do the conversion
                xtemp = orbitalElementConverter(xtemp[0], xtemp[1], xtemp[2], xtemp[3], xtemp[4], xtemp[5]);
                //Save the converted elements
                for(int j=0; j<6; j++) positions[i][j+2] = xtemp[j];
            }
        }
            break;
        default:
            printf("ERROR: Failed to recognize filetype %d\n",filetype);
            break;
    }
    
    //Sanity check section -- print out what we read in
    /*
    cout << "READIN TEST: " << mynames[0] << " " << positions[0][0] << " " << positions[0][1] << " " << positions[0][2] << " " << positions[0][5] << endl;
    cout << "READIN TEST: " << mynames[1] << " " << positions[1][0] << " " << positions[1][1] << " " << positions[1][2] << " " << positions[1][3] << " " << positions[1][4] << " " << positions[1][5] << endl;
     */
    
    return(positions);
}