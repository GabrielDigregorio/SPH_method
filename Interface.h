///**************************************************************************
/// HEADER: Function That Load/Read Input Files And Generate Output Files
///**************************************************************************

#ifndef INTERFACE_H
#define INTERFACE_H
#include "Structures.h"
#include <map>
#include <ctime>
#include <chrono>
#include <cstdlib>
#include <sstream>
#include <iomanip>


// inputReader.cpp
void readParameter(char* parameterFilename, Parameter* parameter);
void readGeometry(char* geometryFilename, Field* field);
void readBrick(int type, std::ifstream* inFile, Field* currentField);

// writeField.cpp
void writeField(Field* field, double t, Format myFormat,
                std::string const &parameterFilename="Undefined",
                std::string const &geometryFilename="Undefined",
                std::string const &filename="result");


// ParaView.cpp
void paraView(std::string const &filename,
              int step,
              std::vector<double> const &pos,
              std::map<std::string, std::vector<double> *> const &scalars,
              std::map<std::string, std::vector<double> *> const &vectors);

void matlab(std::string const &filename,
              std::string const &parameterFilename,
              std::string const &geometryFilename,
              int step,
              std::vector<double> const &pos,
              std::vector<double> const &speed,
              std::vector<double> const &density,
              std::vector<double> const &pressure,
              std::vector<double> const &mass);


#endif
