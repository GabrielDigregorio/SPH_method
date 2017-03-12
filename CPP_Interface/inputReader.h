#ifndef INPUTREADER_H
#define INPUTREADER_H

#include <string>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <map>
#include <cmath>
#include <cassert>
#include <ctime>
#include <fstream>
#include <random>
#include <sstream>
#include <float.h>
#include <iomanip>
// NECESSARY HERE
#include <sstream>
#include <algorithm>

#define N_UL 3
#define N_DATA 9
#define N_PARAM 13

enum geomType{cube,cylinder,sphere};
enum boundCondition{freePart, movingPart, fixedPart};

struct Parameter {
    double kh, k, T, densityRef, B, gamma, g, writeInterval;
    std::string integrationMethod, densityInitMethod, stateEquationMethod, massInitMethod, speedLaw;
};

// Structure qui doit être remplie lors de la lecture du fichier de géométrie (il faudra surement changer la place de cette déclaration aussi). Cette structure contient toute l'information utile de nos simulations.
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
// inputReader.cpp
void readParameter(std::string parameterFilename, Parameter* parameter); ///////////
void readGeometry(std::string geometryFilename, Field* currentField); ////////////
void readBrick(int type, std::ifstream* inFile, Field* currentField, std::vector<int>* sValues); ///////////
// Geometry.cpp
void meshcube(double o[3], double L[3], double s, std::vector<double> &pos, double perturbation = 0.0, bool stack = false);
void meshcylinder(double o[3], double L[3], double s, std::vector<double> &pos, double perturbation = 0.0, bool stack = false);
void meshsphere(double o[3], double L[3], double s, std::vector<double> &pos, double perturbation = 0.0, bool stack = false);
#endif
