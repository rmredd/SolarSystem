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

using std::vector;
using std::valarray;
using std::fstream;
using std::ifstream;
using std::fscanf;
using std::string;

vector<valarray<long double> > ReadParameterFile(string filename, int filetype, vector<string> & names) {
    //Opening up the file
    //Note that '#' is used as a comment character
    ifstream input;
    input.open(filename);
    string myline;
    int nlines = 0;
    while(!input.eof()){
        getline(input,myline);
        if(myline.find("#")==0) continue;
        nlines++;
    }
    input.close();
    
    //Note that the positions valarray contains mass, radius, followed by position information
    vector<string> mynames(nlines);
    vector<valarray<long double> > positions(nlines);
    int place, place_end;
    switch (filetype) {
        case 0:
        {
            //Name(string), X, Y, Z, VX, VY, VZ position data
            input.open(filename);
            valarray<long double> xtemp(6);
            for(int i=0; i<nlines; i++) {
                getline(input,myline);
                //Extract the name of the planet first
                mynames[i] = myline.substr(0,myline.find(" "));
                myline.erase(0,myline.find(" "));
                while(myline.find(" ")==0) myline.erase(0,1);
                place=0;
                for(int j=0; j<6; j++){
                    place_end = myline.find(" ",place);
                    if(place_end==string::npos) place_end = myline.length()-1;
                    xtemp[j] = atof(myline.substr(place,place_end).c_str());
                    place = place_end;
                    while(myline.find(" ",place)==place){
                        place++;
                    }
                }
                positions[i] = xtemp;
            }
            input.close();
        }
            break;
        case 1:
            //Orbital elements -- will need to be converted at some point
        {
            printf("WARNING: Orbital element readin not yet implemented");
        }
            break;
        default:
            printf("ERROR: Failed to recognize filetype %d\n",filetype);
            break;
    }
    
    //Sanity check section -- print out what we read in
    cout << "READIN TEST: " << mynames[0] << " " << positions[0][0] << " " << positions[0][1] << " " << positions[0][2] << " " << positions[0][5] << endl;
    cout << "READIN TEST: " << mynames[1] << " " << positions[1][0] << " " << positions[1][1] << " " << positions[1][2] << " " << positions[1][3] << " " << positions[1][4] << " " << positions[1][5] << endl;
    
    return(positions);
}