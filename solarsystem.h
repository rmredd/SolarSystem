/*
 *  solarsystem.h
 *  SolarSystem
 *
 *  Created by Jadzia on 1/1/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <valarray>

using std::cout;
using std::endl;
using std::vector;
using std::valarray;
using std::string;
using std::ofstream;

#define G 6.67428e-11 //m^3/(kg s^2)
#define Msun 1.989e30 //kg

class Planet {
public:
	void Initialize(double, valarray<double>, double);
	void Iterate(valarray<double>, double);
	double CurrentTime(void) {return(planetTime[planetTime.size()-1]);};
	valarray<double> CurrentPosition(void) {return(X[X.size()-1]);};
	valarray<double> Gravity(const valarray<double>, const double);
	
	void PrintOutput(char*);
	
private:
	double planetMass;  //kg
	vector<double> planetTime;  //time in s
	vector<valarray<double> > X;  //first 3 values x, y, z; next 3 values velocity, m/s
};

void Planet::Initialize(double mass, valarray<double> initialPosition, double initialTime) {
	planetMass = mass;
	X.push_back(initialPosition);
	planetTime.push_back(initialTime);
}

void Planet::Iterate(valarray<double> pos, double time) {
	X.push_back(pos);
	planetTime.push_back(time);
}

valarray<double> Planet::Gravity(const valarray<double> xold, const double t) {
	valarray<double> DX(6);
	
	valarray<double> pPos = X[X.size()-1];
	const double x = xold[0]-pPos[0]; const double y = xold[1]-pPos[1]; const double z = xold[2]-pPos[2];
	const double vx = xold[3]; const double vy = xold[4]; const double vz = xold[5];
	
	//position derivative is velocity
	DX[0] = vx; DX[1] = vy; DX[2] = vz;
	
	//taking velocity derivative
	DX[3] = -G*planetMass/pow( (x*x+y*y+z*z), 1.5)*x;
	DX[4] = -G*planetMass/pow( (x*x+y*y+z*z), 1.5)*y;
	DX[5] = -G*planetMass/pow( (x*x+y*y+z*z), 1.5)*z;
	return(DX);
}

void Planet::PrintOutput(char* filename) {
	int nSteps = X.size();
	ofstream fout;
	fout.open(filename);
	for(int i=0; i<nSteps; i++) {
		fout<<planetTime[i]<<" ";
		fout<<X[i][0]<<" "<<X[i][1]<<" "<<X[i][2]<<" "<<X[i][3]<<" "<<X[i][4]<<" "<<X[i][5]<<" ";
		
		fout<<endl;
	}
	
	fout.close();
}