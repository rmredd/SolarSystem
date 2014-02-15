/*
 *  orbital.h
 *  SolarSystem
 *
 *  Created by Jadzia on 8/5/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "solarsystem.h"

//gives acceleration from Sun's gravitational potential
valarray<double> gravitySun(const valarray<double> xold, const double t) {
	valarray<double> DX(6);
	const double x = xold[0]; const double y = xold[1]; const double z = xold[2];
	const double vx = xold[3]; const double vy = xold[4]; const double vz = xold[5];
	
	//position derivative is velocity
	DX[0] = vx; DX[1] = vy; DX[2] = vz;
	
	//taking velocity derivative
	DX[3] = -G*Msun/pow( (x*x+y*y+z*z), 1.5)*x;
	DX[4] = -G*Msun/pow( (x*x+y*y+z*z), 1.5)*y;
	DX[5] = -G*Msun/pow( (x*x+y*y+z*z), 1.5)*z;
	return(DX);
}


//rotation functions
void rotateZ(double& x, double& y, double& z, double theta) {
	double xold = x; double yold = y;
	x = xold*cos(theta)-yold*sin(theta);
	y = xold*sin(theta)+yold*cos(theta);
}

void rotateX(double& x, double& y, double& z, double theta) {
	double yold = y; double zold = z;
	y = yold*cos(theta)-zold*sin(theta);
	z = yold*sin(theta)+zold*cos(theta);
}

//converts orbital elements and angle to position and velocity
//all input angles in degrees
valarray<double> orbitalElementConverter(double semimajor, double semiminor, double inclination, double ascendingNode, double argumentPerihelion, double angle) {

	valarray<double> myX(6);
	double eccentricity = sqrt(1-semiminor*semiminor/semimajor/semimajor);
	double xSun = eccentricity*semimajor;
	double ySun = 0; double zSun = 0;
	
	double x = 1/sqrt(cos(angle*M_PI/180)*cos(angle*M_PI/180)/semimajor/semimajor+sin(angle*M_PI/180)*sin(angle*M_PI/180)/semiminor/semiminor)*cos(angle*M_PI/180);
	double y = 1/sqrt(cos(angle*M_PI/180)*cos(angle*M_PI/180)/semimajor/semimajor+sin(angle*M_PI/180)*sin(angle*M_PI/180)/semiminor/semiminor)*sin(angle*M_PI/180);
	double z = 0;

	//get position, direction of velocity
	double slope = x/y*semiminor*semiminor/semimajor/semimajor;
	double theta = atan(slope);
	if(y==0) {theta=M_PI/2;}
	double vx = cos(theta);
	double vy = sin(theta);
	if(0 <= angle && angle < 90) {
		vx = -fabs(vx);
		vy = fabs(vy);
	} else if (90<=angle && angle<180) {
		vx = -fabs(vx);
		vy = -fabs(vy);
	} else if (180<=angle && angle<270) {
		vx = fabs(vx);
		vy = -fabs(vy);
	} else {
		vx = fabs(vx);
		vy = fabs(vy);
	}
	double vz = 0;
	
	//calculate energy per mass of satellite
	double EPerMass = G*Msun/semimajor/(1-eccentricity)-1/2*G*Msun/(semimajor-xSun);
	double radius = sqrt((x-xSun)*(x-xSun)+(y-ySun)*(y-ySun));
	double vmag = sqrt(2*(EPerMass+G*Msun/radius));
	cout<<"En v: "<<vmag<<endl;
	cout<<eccentricity<<endl;

	//calculate angular momentum per mass
	double vperi2 = 2*G*Msun/semimajor/(1-eccentricity);
	vperi2 = G*Msun/semimajor*(1+eccentricity)/(1-eccentricity);
	double LPerMass = sqrt(vperi2)*(semimajor-xSun);
	
	//calculate magnitude of velocity
	//mag. of cross product of radius to sun and velocity
	double xR = x-xSun;
	double yR = y-ySun;
	double mag_cross = fabs(vx*yR-vy*xR)/sqrt(xR*xR+yR*yR);
	double mag_v = LPerMass/mag_cross/sqrt(xR*xR+yR*yR);
	cout<<"L v: "<<mag_v<<endl;
	vx *= mag_v;
	vy *= mag_v;
	
	//perform rotations for ascending node and inclination
	//cout<<"test: "<<sqrt(x*x+y*y+z*z)<<endl;
	
	double delta = ascendingNode - argumentPerihelion;
	rotateZ(x,y,z,-delta*M_PI/180); rotateZ(vx,vy,vz,-delta*M_PI/180); rotateZ(xSun,ySun,zSun,-delta*M_PI/180);
	rotateX(x,y,z,inclination*M_PI/180); rotateX(vx,vy,vz,inclination*M_PI/180);
	
	rotateX(xSun,ySun,zSun,inclination*M_PI/180);
	rotateZ(x,y,z,delta*M_PI/180); rotateZ(vx,vy,vz,delta*M_PI/180); rotateZ(xSun,ySun,zSun,delta*M_PI/180);
	//cout<<"test: "<<sqrt(x*x+y*y+z*z)<<endl;

	
	//perform rotation for argument of perihelion
	rotateZ(x,y,z,argumentPerihelion*M_PI/180); rotateZ(vx,vy,vz,argumentPerihelion*M_PI/180);
	rotateZ(xSun,ySun,zSun,argumentPerihelion*M_PI/180);
	//cout<<"test: "<<sqrt(x*x+y*y+z*z)<<endl;

	//set Sun to origin and set values
	myX[0] = x-xSun; myX[1] = y-ySun; myX[2] = z-zSun; myX[3] = vx; myX[4] = vy; myX[5] = vz;
	
	return(myX);
}