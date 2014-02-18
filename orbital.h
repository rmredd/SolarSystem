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
valarray<long double> gravitySun(const valarray<long double> xold, const long double t) {
	valarray<long double> DX(6);
	const long double x = xold[0]; const long double y = xold[1]; const long double z = xold[2];
	const long double vx = xold[3]; const long double vy = xold[4]; const long double vz = xold[5];
	
	//position derivative is velocity
	DX[0] = vx; DX[1] = vy; DX[2] = vz;
	
	//taking velocity derivative
	DX[3] = -G*Msun/pow( (x*x+y*y+z*z), 1.5)*x;
	DX[4] = -G*Msun/pow( (x*x+y*y+z*z), 1.5)*y;
	DX[5] = -G*Msun/pow( (x*x+y*y+z*z), 1.5)*z;
	return(DX);
}


//rotation functions
void rotateZ(long double& x, long double& y, long double& z, long double theta) {
	long double xold = x; long double yold = y;
	x = xold*cos(theta)-yold*sin(theta);
	y = xold*sin(theta)+yold*cos(theta);
}

void rotateX(long double& x, long double& y, long double& z, long double theta) {
	long double yold = y; long double zold = z;
	y = yold*cos(theta)-zold*sin(theta);
	z = yold*sin(theta)+zold*cos(theta);
}

//converts orbital elements and angle to position and velocity
//all input angles in degrees
valarray<long double> orbitalElementConverter(long double semimajor, long double semiminor, long double inclination, long double ascendingNode, long double argumentPerihelion, long double angle) {

	valarray<long double> myX(6);
	long double eccentricity = sqrt(1-semiminor*semiminor/semimajor/semimajor);
	long double xSun = eccentricity*semimajor;
	long double ySun = 0; long double zSun = 0;
	
	long double x = 1/sqrt(cos(angle*M_PI/180)*cos(angle*M_PI/180)/semimajor/semimajor+sin(angle*M_PI/180)*sin(angle*M_PI/180)/semiminor/semiminor)*cos(angle*M_PI/180);
	long double y = 1/sqrt(cos(angle*M_PI/180)*cos(angle*M_PI/180)/semimajor/semimajor+sin(angle*M_PI/180)*sin(angle*M_PI/180)/semiminor/semiminor)*sin(angle*M_PI/180);
	long double z = 0;

	//get position, direction of velocity
	long double slope = x/y*semiminor*semiminor/semimajor/semimajor;
	long double theta = atan(slope);
	if(y==0) {theta=M_PI/2;}
	long double vx = cos(theta);
	long double vy = sin(theta);
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
	long double vz = 0;
	
	//calculate energy per mass of satellite
	long double EPerMass = G*Msun/semimajor/(1-eccentricity)-1/2*G*Msun/(semimajor-xSun);
	long double radius = sqrt((x-xSun)*(x-xSun)+(y-ySun)*(y-ySun));
	long double vmag = sqrt(2*(EPerMass+G*Msun/radius));
	cout<<"En v: "<<vmag<<endl;
	cout<<eccentricity<<endl;

	//calculate angular momentum per mass
	long double vperi2 = 2*G*Msun/semimajor/(1-eccentricity);
	vperi2 = G*Msun/semimajor*(1+eccentricity)/(1-eccentricity);
	long double LPerMass = sqrt(vperi2)*(semimajor-xSun);
	
	//calculate magnitude of velocity
	//mag. of cross product of radius to sun and velocity
	long double xR = x-xSun;
	long double yR = y-ySun;
	long double mag_cross = fabs(vx*yR-vy*xR)/sqrt(xR*xR+yR*yR);
	long double mag_v = LPerMass/mag_cross/sqrt(xR*xR+yR*yR);
	cout<<"L v: "<<mag_v<<endl;
	vx *= mag_v;
	vy *= mag_v;
	
	//perform rotations for ascending node and inclination
	//cout<<"test: "<<sqrt(x*x+y*y+z*z)<<endl;
	
	long double delta = ascendingNode - argumentPerihelion;
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