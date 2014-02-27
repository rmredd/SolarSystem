//
//  parameter_readin.h
//  SolarSystem
//
//  Created by Rachel Reddick on 2/22/14.
//
//

//General function that manages all our basic readin stuffs

#include valarray
#include vector
#include ifstream

using std::vector;
using std::valarray;
using std::ifstream;

vector<valarray<long double> > ReadParameterFile(string filename, int filetype, vector<string> & names) {
    //Opening up the file
    ifstream input (filename);
    string myline;
    int nlines = 0;
    while(!input.eof()){
        myline = input.getline();
        nlines++;
    }
    input.close();
    
    //Note that the positions valarray contains mass, radius, followed by position information
    
    vector<valarray<long double> > positions;
    
    return(positions);
}