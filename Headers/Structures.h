///**************************************************************************
/// HEADER: Contains All Structures
///**************************************************************************

#ifndef STRUCTURES_H
#define STRUCTURES_H
#include "Main.h"

// Error types
enum Error {noError, argumentError, parameterError, geometryError, consistencyError, NB_ERROR_VALUE};
// Kernel Types
enum Kernel {Gaussian, Bell_shaped, Cubic_spline, Quadratic, Quintic, Quintic_spline, NB_KERNEL_VALUE};

// model of viscosity formulation
enum ViscosityModel {violeauArtificial, NB_VISCOSITY_VALUE};

// integrationMethod = euler ou RK2
enum IntegrationMethod {euler, RK2, NB_INTEGRATION_VALUE};

// AdaptativeTimeStep
enum AdaptativeTimeStep {no, yes, NB_ADAPTATIVE_VALUE};

// densityInitMethod = hydrosatic, etc.
enum DensityInitMethod {hydrostatic, homogeneous, NB_DENSITYINIT_VALUE};

// stateEquationMethod = quasiIncompressible, perfectGas, etc.
enum StateEquationMethod {quasiIncompressible, perfectGas, NB_STATEEQUATION_VALUE};

// massInitMethod = violeau2012 (all particles have same volumes), etc.
enum MassInitMethod {violeau2012, NB_MASSINIT_VALUE};

// speedLaw = Will dictate the behaviour of moving boundaries: constant, sine, exponential
enum SpeedLaw {constant, sine, exponential,level_arm, NB_SPEEDLAW_VALUE};

// Write Format output
enum Format {ParaView, Matlab, Both, NB_FORMAT_VALUE};

// charactTime = characteristic time of movement or period of oscillations

// movingDirection = direction de la paroi mouvante

/* Parameter Structure
 * kh = smothing length
 * k = time step
 * T = simulation time
 * densityRef = density of the fluid at atmospheric pressure
 * l & u = lower and upper limit of the domain
 * B & gamma = fluid constants
 * g = gravity
 * writeInteval = time interval between two outputs file are generated
 * c, alpha, beta, epsilon = parameter relativ to artificial viscosity
*/


struct Parameter {
    double kh, k, T, densityRef, B, gamma, g, writeInterval, charactTime, c, alpha, beta, epsilon, molarMass, temperature, theta;
    double movingDirection[3];
    Kernel kernel;
    ViscosityModel viscosityModel;
    IntegrationMethod integrationMethod;
    AdaptativeTimeStep adaptativeTimeStep;
    DensityInitMethod densityInitMethod;
    StateEquationMethod stateEquationMethod;
    MassInitMethod massInitMethod;
    SpeedLaw speedLaw;
    Format format;
};


struct Field {
    int nFree, nFixed, nMoving, nTotal;

    double l[3];
    double u[3];

    double nextK=0.0;

    std::vector<double> pos;

    std::vector<double> speed;

    std::vector<double> density;

    std::vector<double> pressure;

    std::vector<double> mass;

    std::vector<double> info_block;

    std::vector<double> info_moving;

};



#endif
