#include "Main.h"
#include "Interface.h"

/*struct Parameter {
    double kh, k, T, densityRef, B, gamma, g, writeInterval, charactTime, c, alpha, beta, epsilon;
    double movingDirection[3];
    Kernel kernel;
    ViscosityModel viscosityModel;
    IntegrationMethod integrationMethod;
    DensityInitMethod densityInitMethod;
    StateEquationMethod stateEquationMethod;
    MassInitMethod massInitMethod;
    SpeedLaw speedLaw;
    Format format;
    }
*/

void Consistency(Parameter* param,Field* field)
{
double ux=field->u[0];
double uy=field->u[1];
double uz=field->u[2];
double lx=field->l[0];
double ly=field->l[1];
double lz=field->l[2];
int nFree=field->nFree;
int nFixed=field->nFixed;
int nMoving= field->nMoving;
int nTotal=field->nTotal;
// // check of structure "Field"
// check domain (l vs u)
assert(ux>lx && uy>ly && uz>lz && "ui not > li");
// check number of particules
assert(nFree>=0 && nFixed>=0 && nMoving>=0 && nTotal>=0 && "number of particule not an positive integer ");
// check position of the cube, boundary... if include in the domain ??

// check s for free, moving fixed.



double kh=param->kh;
double k=param->k;
double T=param->T;
double densityRef=param->densityRef;
double B=param->B;
double gamma=param->gamma;
double writeInterval=param->writeInterval;
double alpha=param->alpha; 
double beta=param->beta;
double epsilon=param->epsilon;
double 
// // check of struture "Parameter"
// check kh >0 (and kh>1.2*s ??)
assert(kh>0.0 );
// check k>0 && k< T
assert(k>0.0 && k<T && "time step not good " );
// check density initial
assert(densityRef>0.0 && "density not positive");
// check B >0 , and what else ?
assert(B>0.0 && "B is not positive");
// check gamma >0 , and what else ?
assert(gamma>0 && "gamma is not positive");
// check write interval >k and <T
assert(writeInterval>k && writeInterval<T && "writeInterval not good");
// check  charactTime is positive 
assert(charactTime>=0 && "charactTime is not positive");
// check c ( speed of sound ) is positive
assert(c>0 && "c is not positive");
// check alpha is a positve number
assert(alpha>=0 && "alpha is not positive");
// check beta is a positve number
assert(beta>=0 && "beta is not positive");
// check epsilon is a positive number;
assert(epsilon>0 && "beta is not positive");



// // check Consistency

}