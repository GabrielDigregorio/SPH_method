///**************************************************************************
/// HEADER: Contains All Structures
///**************************************************************************

#ifndef STRUCTURES_H
#define STRUCTURES_H


// Kernel Types
enum Kernel { Gaussian=1, Bell_shaped=2, Cubic_spline=3, Quadratic=4, Quintic=5, Quintic_spline=6 };



/* Parameter Structure
 * kh = smothing length
 * k = time step
 * T = simulation time
 * densityRef = density of the fluid at atmospheric pressure
 * l & u = lower and upper limit of the domain
 * B & gamma = fluid constants
 * g = gravity
 * writeInteval = time interval between two outputs file are generated
 * integrationMethod = euler ou RK2
 * densityInitMethod = hydrosatic, etc.
 * stateEquationMethod = quasiIncompressible, perfectGas, etc.
 * massInitMethod = violeau2012 (all particles have same volumes), etc.
 * speedLaw = To be determined, will dictate the behaviour of moving boundaries
*/

struct Parameter {
    double kh, k, T, densityRef, B, gamma, g, writeInterval;
    std::string integrationMethod, densityInitMethod, stateEquationMethod, massInitMethod, speedLaw;
};

struct Field {
    std::vector<double> sFree;
    std::vector<double> sMoving;
    std::vector<double> sFixed;

    double l[3];
    double u[3];

    std::vector<double> posFree;
    std::vector<double> posMoving;
    std::vector<double> posFixed;

    std::vector<double> speedFree;
    std::vector<double> speedMoving;
    //Speed fixed = 0 of course

    std::vector<double> densityFree;
    std::vector<double> densityMoving;
    std::vector<double> densityFixed;

    std::vector<double> pressureFree;
    std::vector<double> pressureMoving;
    std::vector<double> pressureFixed;

    std::vector<double> massFree;
    std::vector<double> massMoving;
    std::vector<double> massFixed;
};

#endif
