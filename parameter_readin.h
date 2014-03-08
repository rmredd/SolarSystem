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
    vector<int> is_comment;
    int nlines = 0; int nplanets = 0;
    while(!input.eof()) {
        getline(input,myline);
        nlines++;
        if(myline.find("#")==0) {
            is_comment.push_back(1);
        } else {
            nplanets++;
            is_comment.push_back(0);
        }
    }
    input.close();
    
    vector<valarray<long double> > positions(nplanets);
    valarray<long double> xtemp(nentries);
    vector<string> mynames(nplanets);
    
    cout << "Reading in " << nplanets << " objects" << endl;
    //for(int i=0; i<nlines; i++) cout << is_comment[i] << " ";
    //cout << endl;
    
    int i_planet = 0, place, place_end, name_end;
    input.open(filename);
    for(int i=0; i<nlines; i++) {
        //cout << "LINE READ: " << i << " of " << nlines << " " << is_comment[i] << endl;
        getline(input,myline);
        //Check that this isn't a comment
        if(is_comment[i]==1) {
            //cout << "COMMENT!!" << endl;
            continue;
        }
        
        //Extract the name of the planet first
        //for(int j=0; j<nlines; j++) cout << is_comment[j] << " ";
        //cout << endl;
        name_end = myline.find(" ");
        mynames[i_planet] = myline.substr(0,name_end);
        //for(int j=0; j<nlines; j++) cout << is_comment[j] << " ";
        //cout << endl;
        cout << mynames[i_planet] << endl;
        myline.erase(0,name_end);
        //cout << myline << endl;
        while(myline.find(" ")==0) {
            myline.erase(0,1);
        }
        place=0;
        for(int j=0; j<nentries; j++){
            //cout << "    Entry: " << j << " of " << nentries << endl;
            place_end = myline.find(" ",place);
            if(place_end==string::npos) place_end = myline.length()-1;
            if(j<7 && place_end == myline.length()-1) {
                cout << "WARNING: Error in file readin at line " << i << endl;
                cout << "         Insufficient entries in line for chosen format" << endl;
                break;
            }
            xtemp[j] = atof(myline.substr(place,place_end).c_str());
            place = place_end;
            //cout << "    " << j << " " << xtemp[j] << endl;
            while(myline.find(" ",place)==place){
                place++;
            }
        }
        //for(int j=0; j<nentries; j++) cout << xtemp[j] << " ";
        //cout << endl;
        positions[i_planet] = xtemp;
        //cout << i << " " << xtemp[0] << endl;
        i_planet++;
        //cout << "Out of the loop" << endl;
    }
    cout << positions[0][0] << " " << positions[1][0] << endl;
    cout << "Finished reading" << endl;
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