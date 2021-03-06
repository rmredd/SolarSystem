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
#include <string>

using std::cout;
using std::endl;
using std::vector;
using std::valarray;
using std::string;
using std::ofstream;
using std::string;

//Including the integrators

#define G 6.67384e-11 //m^3/(kg s^2)
#define Msun 1.989e30 //kg

class Planet {
public:
	void Initialize(long double, long double, valarray<long double>, long double, string);
	void Iterate(valarray<long double>, valarray<long double>, long double);
	long double CurrentTime(void) {return(planetTime);};
    long double Radius(void) {return(planetRadius);};
    long double Mass(void) {return(planetMass);};
    valarray<long double> CurrentPosition(void) {return(X);};
    valarray<long double> CurrentPosition(int); //Returns position at specified step
    valarray<long double> CurrentDerivative(int); //Returns derivative at specified step, if already evaluated
    void SaveFuturePosition(valarray<long double>); //Save a future position
	valarray<long double> Gravity(const valarray<long double>, const long double); //Gravity acting on another body
    valarray<long double> PositionDerivative(void); //Velocity derivative part
    
    string MyNameIs(void) {return(Name);} //returns planet's name as a string
	
private:
	long double planetMass;  //kg
	long double planetTime;  //current time in s
    long double planetRadius; //planet's radius in m
	valarray<long double> X;  //first 3 values x, y, z; next 3 values velocity, m/s
    valarray<long double> X1; //location at current timestep-1
    valarray<long double> X2; //location at current timestep-2
    valarray<long double> X3; //location at current timestep-3
    valarray<long double> Xfuture; //location at future timestep+1
    //Derivatives previously evaluated
    valarray<long double> DX;
    valarray<long double> DX1;
    valarray<long double> DX2;
    valarray<long double> DX3;
    string Name; //Name of the planet
};

//Initialize your planet!
void Planet::Initialize(long double mass, long double radius, valarray<long double> initialPosition, long double initialTime, string myName) {
	planetMass = mass;
    planetRadius = radius;
	X = initialPosition;
	planetTime = initialTime;
    X1 = initialPosition;
    X2 = initialPosition;
    X3 = initialPosition;
    Xfuture = initialPosition;
    Name = myName;
}

void Planet::Iterate(valarray<long double> pos, valarray<long double> der, long double time) {
    X3 = X2;
    X2 = X1;
    X1 = X;
	X = pos;
    DX3 = DX2;
    DX2 = DX1;
    DX1 = der;
	planetTime = time;
}

//An extra future position that we can save
void Planet::SaveFuturePosition(valarray<long double> xpos){
    Xfuture = xpos;
}

//Returns planet's position.  0 = current step, 1 = one step previous, etc.
//Defaults to current timestep
valarray<long double> Planet::CurrentPosition(int mystep){
    switch(mystep){
        case 0:
            return(X);
        case 1:
            return(X1);
        case 2:
            return(X2);
        case 3:
            return(X3);
        case -1:
            return(Xfuture); //returns a position from the future
        default:
            return(X);
    }
}

//Returns planet's derivative at a given timestep
valarray<long double> Planet::CurrentDerivative(int mystep) {
    switch(mystep){
        case 0:
            return(DX);
        case 1:
            return(DX1);
        case 2:
            return(DX2);
        case 3:
            return(DX3);
        default:
            return(DX);
    }
}

//Effect of a planet's gravity on another object at position x_else
valarray<long double> Planet::Gravity(const valarray<long double> x_else, const long double t) {
	valarray<long double> DX(6);
	
	const long double x = x_else[0]-X[0]; const long double y = x_else[1]-X[1]; const long double z = x_else[2]-X[2];
	//const long double vx = x_else[3]; const long double vy = x_else[4]; const long double vz = x_else[5];
	
	//position derivative is velocity, which is not given by this object
	DX[0] = 0.; DX[1] = 0.; DX[2] = 0.;
	
	//taking velocity derivative
	DX[3] = -G*planetMass/pow( (x*x+y*y+z*z), 1.5)*x;
	DX[4] = -G*planetMass/pow( (x*x+y*y+z*z), 1.5)*y;
	DX[5] = -G*planetMass/pow( (x*x+y*y+z*z), 1.5)*z;
	return(DX);
}

//Derivative -- note that this is constant velocities in the absence of other objects
valarray<long double> Planet::PositionDerivative(){
    valarray<long double> DX(6);
    DX[0] = X[3]; DX[1] = X[4]; DX[2] = X[5];
    DX[3] = 0.; DX[4] = 0.; DX[5] = 0.;
    return(DX);
}

//Main class for operating on everything
class SolarSystem {
public:
    void Initialize(long double, string); //Initialize -- all this does is set the timestep size and outputs
    void AddPlanet(long double, long double, valarray<long double>, long double, string); //Routine for adding a planet -- repeat to fully initialize
    void PrintPlanets(void); //Append the position of the current timestep to files
    int CheckForCollision(bool); //Test to see if worlds collided

    vector<valarray<long double> > Gravity(int); //Calculates dervative at current time step
    valarray<long double> Gravity(const valarray<long double>, const long double); //Same, but works with Runge Kutta
    
    //Utility function for centering on the system's center of mass
    void MoveToCenterOfMass(void);
    
    //Some input-and-output array creation items necessary for iteration setup
    int NumberOfPlanets(void) {return(nPlanets);};
    int NumberOfSteps(void) {return(nSteps);};
    long double MyTimestep(void) {return(Timestep);};
    valarray<long double> CurrentPosition(int); //Returns current positions for one planet
    vector<valarray<long double> > CurrentPositionAllPlanets(int); //Returns indicated timestep position for all planets
    void UpdateFuturePositions(vector<valarray<long double> >); //Updates future positions for all planets
    void Iterate(vector<valarray<long double> >, vector<valarray<long double> >); //Iterates, taking in the new position as current and previous step's derivative
    
private:
    vector<Planet> Planets;
    long double Timestep; //timestep
    string Output_dir; //output directory
    vector<string> OutputFiles; //Files to be used for output, one per planet
    int nPlanets; //a count of how many planets we've got
    long nSteps; //Number of timesteps we've taken so far
};

//Basic initialization step
void SolarSystem::Initialize(long double dt, string outdir){
    Timestep = dt;
    Output_dir = outdir;
    nPlanets = 0;
    nSteps = 0;
}

//Make a planet, initialize it, and add it to our list of planets
//Note that the input is assumed to be an X, Y, Z, VX, VY, VZ vector
void SolarSystem::AddPlanet(long double mass, long double radius, valarray<long double> initialPosition, long double initialTime, string myName){
    
    Planet newPlanet;
    newPlanet.Initialize(mass,radius,initialPosition, initialTime, myName);
    Planets.push_back(newPlanet);
    OutputFiles.push_back(Output_dir+myName+".dat");
    nPlanets++;
}

//Function for appending current positions and velocities of planets to separate output files,
//Which are determined by their names
void SolarSystem::PrintPlanets(){
    ofstream fs;
    valarray<long double> position(6);
    for(int i=0; i<nPlanets; i++){
        if(nSteps>0) {
            fs.open(OutputFiles[i],ofstream::app);
        } else {
            fs.open(OutputFiles[i]);
        }
        fs.precision(10);
        position = Planets[i].CurrentPosition(0);
        fs << position[0] << " " << position[1] << " " <<position[2] << " " << position[3] << " " << position[4] << " " << position[5] << endl;
        fs.close();
    }
    //And fill in a separate file with time information
    if(nSteps > 0) {
        fs.open(Output_dir+"time.dat",ofstream::app);
    } else {
        fs.open(Output_dir+"time.dat");
    }
    double mytime;
    mytime = Planets[0].CurrentTime();
    fs << nSteps << " " << mytime << " " << mytime/3600./24./365.24 << endl;
}

//Function that tests for collisions between planets
//Prints an alert if there is a collision and returns true; returns false otherwise
int SolarSystem::CheckForCollision(bool halt_for_projected_collision){
    long double separation, sum_radii;
    valarray<long double> p1(6), p2(6);
    for(int i=0; i<nPlanets; i++){
        p1 = Planets[i].CurrentPosition(0);
        for(int j=i+1; j<nPlanets; j++){
            p2 = Planets[j].CurrentPosition(0);
            //Calculate the distance between the two planets
            separation = sqrt( (p1[0]-p2[0])*(p1[0]-p2[0]) + (p1[1]-p2[1])*(p1[1]-p2[1]) + (p1[2]-p2[2])*(p1[2]-p2[2]));
            //Is this distance less than the sum of their radii?
            sum_radii = Planets[i].Radius() + Planets[j].Radius();
            if (separation < sum_radii) {
                //There is a collision!  Print a warning and return true
                cout << "COLLISION ALERT!  Planets " << Planets[i].MyNameIs() << " and " << Planets[j].MyNameIs() << " collide at time " << Planets[i].CurrentTime() << endl;
                return(1);
            }
            if (separation > sum_radii && separation < 20*sum_radii && halt_for_projected_collision) {
                //There is a potential collision here -- let's do a more complex check
                long double close_approach, E, Lmag2, Gmm;
                valarray<long double> CoM(6), L(3);
                //Get the center of mass positon and motion, and subtract it off
                for(int k=0; k<6; k++) {
                    CoM[k] = (Planets[i].Mass()*p1[k]+Planets[j].Mass()*p2[k])/(Planets[i].Mass() + Planets[j].Mass());
                }
                p1 = p1 - CoM;
                p2 = p2 - CoM;
                //Get the orbital energy
                Gmm = G*Planets[i].Mass()*Planets[j].Mass();
                E = 0.5*Planets[i].Mass()*(p1[3]*p1[3]+p1[4]*p1[4]+p1[5]*p1[5]) + 0.5*Planets[j].Mass()*(p2[3]*p2[3]+p2[4]*p2[4]+p2[5]*p2[5]) - Gmm/separation;
                //Calculating the angular momentum vector
                for(int k=0; k<3; k++) L[k] = Planets[i].Mass()*(p1[(1+k)%3]*p1[3+((2+k)%3)] - p1[(2+k)%3]*p1[3+(1+k)%3]) + Planets[j].Mass()*(p2[(1+k)%3]*p2[3+((2+k)%3)] - p2[(2+k)%3]*p2[3+(1+k)%3]);
                //Now get the amplitude of the vector
                Lmag2 = L[0]*L[0]+L[1]*L[1]+L[2]*L[2];
                //Now calculate the separation at closest approach
                close_approach = (-Gmm + sqrt(Gmm*Gmm + 2*E*Lmag2*(1/Planets[i].Mass() + 1/Planets[j].Mass())))/(2*E);
                if(close_approach < sum_radii) {
                    cout << "COLLISION WARNING!  Planets " << Planets[i].MyNameIs() << " and " << Planets[j].MyNameIs() << " are on a collision course at time " << Planets[i].CurrentTime() << endl;
                    cout << close_approach << " " << sum_radii << " " << E << " " << Lmag2 << endl;
                    return(2);
                }
            }
        }
    }
    return(0);
}

//Make the Gravity function, which returns the derivative we need at the current time step
//Note this works at discrete steps, rather than doing the operation between, as per
//the Runge Kutta integrator
vector<valarray<long double> > SolarSystem::Gravity(int mystep){
    vector<valarray<long double> > DX(nPlanets);

    //Initialize our arrays
    valarray<long double> temp(6), temp2(6);
    long double x, y, z;

    //Case where we don't have the derivative calculated yet
    if(mystep<=0) {
        for(int i=0; i<6; i++) temp[i] = 0;
        for(int i=0; i<nPlanets; i++) DX[i] = temp;
    
        for(int i=0; i<nPlanets; i++){
            for(int j=i; j<nPlanets; j++) {
                if (i==j) {
                    //This is matching the planet against itself; move the velocities to get dx/dt
                    DX[i][0] = Planets[i].CurrentPosition(mystep)[3];
                    DX[i][1] = Planets[i].CurrentPosition(mystep)[4];
                    DX[i][2] = Planets[i].CurrentPosition(mystep)[5];
                } else {
                    //Otherwise, we're comparing two planets.  Update dv/dt for both with gravity
                    temp = Planets[i].CurrentPosition(mystep);
                    temp2 = Planets[j].CurrentPosition(mystep);
                    x = temp[0] - temp2[0];
                    y = temp[1] - temp2[1];
                    z = temp[2] - temp2[2];
                
                    temp[0] = -G/pow( (x*x+y*y+z*z), 1.5)*x;
                    temp[1] = -G/pow( (x*x+y*y+z*z), 1.5)*y;
                    temp[2] = -G/pow( (x*x+y*y+z*z), 1.5)*z;
                    DX[i][3] += Planets[j].Mass()*temp[0];
                    DX[i][4] += Planets[j].Mass()*temp[1];
                    DX[i][5] += Planets[j].Mass()*temp[2];
                    DX[j][3] -= Planets[i].Mass()*temp[0];
                    DX[j][4] -= Planets[i].Mass()*temp[1];
                    DX[j][5] -= Planets[i].Mass()*temp[2];
                }
            }
        }
    } else {
        //We've already computed the derivative, so just pull it from the planets' vectors
        for(int i=0; i<nPlanets; i++) {
            DX[i] = Planets[i].CurrentDerivative(mystep);
        }
    }
    
    //cout << "Gravity test: " << mystep << " " << x << " " << y << " " << z << endl;
    
    return(DX);
}

//A version of the Gravity function that works with the general Runge Kutta integrator
//Note that this still needs the class data for the masses, number of planets, and requires
//Conversion of current position into a single valarray
valarray<long double> SolarSystem::Gravity(const valarray<long double> position, const long double t){
    
    valarray<long double> DX(nPlanets*6);
    valarray<long double> temp(6);
    long double x, y, z;
    for(int i=0; i<6; i++) temp[i] = 0;
    for(int i=0; i<nPlanets*6; i++) DX[i] = 0;
    
    //Run through everything, including the gravity calculation
    for(int i=0; i<nPlanets; i++){
        for(int j=i; j<nPlanets; j++){
            if(i==j){
                //Transferring over the velocities
                DX[6*i] = position[6*i+3];
                DX[6*i+1] = position[6*i+4];
                DX[6*i+2] = position[6*i+5];
            } else {
                //Now getting the accelerations
                x = position[6*i]-position[6*j];
                y = position[6*i+1]-position[6*j+1];
                z = position[6*i+2]-position[6*j+2];
                temp[0] = -G/pow( (x*x+y*y+z*z), 1.5)*x;
                temp[1] = -G/pow( (x*x+y*y+z*z), 1.5)*y;
                temp[2] = -G/pow( (x*x+y*y+z*z), 1.5)*z;
                DX[6*i+3] += Planets[j].Mass()*temp[0];
                DX[6*i+4] += Planets[j].Mass()*temp[1];
                DX[6*i+5] += Planets[j].Mass()*temp[2];
                DX[6*j+3] -= Planets[i].Mass()*temp[0];
                DX[6*j+4] -= Planets[i].Mass()*temp[1];
                DX[6*j+5] -= Planets[i].Mass()*temp[2];
            }
        }
    }
    
    return(DX);
}

//Function for recentering the system, so it doesn't wander too much
//Sets [0, 0, 0] to be stationary center of the system
void SolarSystem::MoveToCenterOfMass() {
    //Set up some variables
    long double total_mass=0;
    valarray<long double> centerOfMass(6);
    for(int i=0; i<6; i++) centerOfMass[i] = 0;
    
    //First, calculate the center of mass of the system
    for(int i=0; i<nPlanets; i++){
        total_mass += Planets[i].Mass();
        centerOfMass += Planets[i].CurrentPosition()*Planets[i].Mass();
    }
    centerOfMass /= total_mass;
    
    //Now, substract the center of mass position and velocity for each planet
    for(int i=0; i<nPlanets; i++){
        Planets[i].Initialize(Planets[i].Mass(), Planets[i].Radius(), Planets[i].CurrentPosition()-centerOfMass, Planets[i].CurrentTime(), Planets[i].MyNameIs());
    }
    
}

//Function for having a look at any one planet's current position
valarray<long double> SolarSystem::CurrentPosition(int myplanet){
    if ( myplanet >= nPlanets || myplanet < 0){
        //Default -- return the position of the zero index planet
        return(Planets[0].CurrentPosition());
    }
    return(Planets[myplanet].CurrentPosition());
}

//Function for having a look at all planets current positions at specified timestep
vector<valarray<long double> > SolarSystem::CurrentPositionAllPlanets(int mystep){
    vector<valarray<long double> > positions(nPlanets);
    for(int i=0; i<nPlanets; i++) {
        positions[i] = Planets[i].CurrentPosition(mystep);
    }
    return(positions);
}

//Function that updates the future positions of all planets
void SolarSystem::UpdateFuturePositions(vector<valarray<long double> > xfuture) {
    for(int i=0; i<nPlanets; i++) {
        Planets[i].SaveFuturePosition(xfuture[i]);
    }
}

//
void SolarSystem::Iterate(vector<valarray<long double> > xnew, vector<valarray<long double> > der_new) {
    nSteps++;
    for(int i=0; i<nPlanets; i++) Planets[i].Iterate(xnew[i], der_new[i], nSteps*Timestep);
}

//Version of the Runge Kutta integrator which operates on a SolarSystem object
//Runge Kutta integration step
valarray<long double> RungeKuttaSystem(SolarSystem &system, valarray<long double> &der) {
    int number_of_planets = system.NumberOfPlanets();
    valarray<long double> temp(6);
    int size = 6*number_of_planets;
    valarray<long double> xold(size), xnew(size);
    //Setting up the input vector
    for(int i=0; i<number_of_planets; i++) {
        temp = system.CurrentPosition(i);
        for(int j=0; j<6; j++) xold[6*i+j] = temp[j];
    }
    valarray<long double> k1(size), k2(size), k3(size), k4(size);
    long double h = system.MyTimestep();
    long double t = h*system.NumberOfSteps();
    der = system.Gravity(xold,t);
    k1=h*der;
    k2=h*system.Gravity(xold+k1,t+h/2);
    k3=h*system.Gravity(xold+(.5)*k2,t+h/2);
    k4=h*system.Gravity(xold+(.5)*k3,t+h);
    xnew = xold+(k1+(2.0)*k2+(2.0)*k3+k4)/6;
    return(xnew);
}

//Version of the predictor-corrector integrator that operates on a SolarSystem object
vector<valarray<long double> > PredictorCorrectorSystem(SolarSystem &system, vector<valarray<long double> > &der){
    int number_of_planets = system.NumberOfPlanets();
    vector<valarray<long double> > xold(number_of_planets), xnew(number_of_planets);
    vector<valarray<long double> > der_xm1(number_of_planets), der_xm2(number_of_planets), der_xm3(number_of_planets), der_xold(number_of_planets), der_xp1(number_of_planets), xp1(number_of_planets);
    
    xold = system.CurrentPositionAllPlanets(0);
    

    der_xold = system.Gravity(0);
    der_xm1 = system.Gravity(1);
    der_xm2 = system.Gravity(2);
    der_xm3 = system.Gravity(3);
    
    //Saving the newest derivative for later use
    der = der_xold;
    
    //Getting the time step
    long double h = system.MyTimestep();
    
    //Adams-Bashford Prediction step
    for(int i=0; i<number_of_planets; i++) xp1[i] = xold[i] + h*(55./24.*der_xold[i] - 59./24.*der_xm1[i] + 37./24.*der_xm2[i] - 0.375*der_xm3[i]);
    
    //Add this future step to the system
    system.UpdateFuturePositions(xp1);
    
    //Evaluation step
    der_xp1 = system.Gravity(-1);
    
    //Adams-Moulton Correction step
    for(int i=0; i<number_of_planets; i++) xnew[i] = xold[i] + h*( 0.375*der_xp1[i] + 19./24.*der_xold[i] - 5./24.*der_xm1[i] + 1./24.*der_xm2[i]);
    
    return(xnew);
}

//Take a deep breath -- we're taking a single time step.  Should do only N(N-1)/2 calcuations
void StepTheSystem(SolarSystem &system) {
    //What we do depends on how many steps have already been taken;
    //if we're on the third step or further along, we can use the PEC algorithm; otherwise,
    //we set up for the PEC using Runge Kutta
    //Here's the vector of valarrays we'll be using to update
    int number_of_planets = system.NumberOfPlanets();
    int nsteps = system.NumberOfSteps();
    
    vector<valarray<long double> > xnew(number_of_planets);
    vector<valarray<long double> > der_new(number_of_planets);
    valarray<long double> temp(6), der_temp(6);
    valarray<long double> end_position(6*number_of_planets), end_der(6*number_of_planets); //Our output vector
    
    if(nsteps < 3){
        //Using Runge-Kutta
        end_position = RungeKuttaSystem(system, end_der);
        
        //And do a quick wrap-around to separate out the different planets
        for(int i=0; i<number_of_planets; i++) {
            for(int j=0; j<6; j++) {
                temp[j] = end_position[6*i+j];
                der_temp[j] = end_der[6*i+j];
            }
            xnew[i] = temp;
            der_new[i] = der_temp;
        }

    
    } else {
        //Using PEC (predictor-corrector)
        xnew = PredictorCorrectorSystem(system, der_new);
        if (nsteps < 10) {
            cout << nsteps << " " << xnew[0][0] << " " << xnew[0][1] << " " << xnew[0][2] << endl;
        }
    }
    
    //Run each planet's internal iterations in order to update correctly
    system.Iterate(xnew, der_new);
}

//Run the whole shebang for a specified period of time (in seconds)
//This includes collision checking and the writing of output for each planet's position
void RunTheSystem(SolarSystem &system, long double max_time, int steps_between_prints, bool halt_for_projected_collision) {
    int collided = 0;
    long nsteps=1;
    while(system.MyTimestep()*system.NumberOfSteps() < max_time && collided==0) {
        StepTheSystem(system);
        //cout << "Running at step: " << system.NumberOfSteps() << endl;
        if(nsteps % steps_between_prints == 0) system.PrintPlanets();
        collided = system.CheckForCollision(halt_for_projected_collision);
        nsteps++;
    }
    //print data for the last step if there's a collision
    if(collided>0) system.PrintPlanets();
}

