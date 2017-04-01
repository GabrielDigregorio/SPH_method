///**************************************************************************
/// HEADER: Function That Load/Read Input Files And Generate Output Files
///**************************************************************************

#ifndef INTERFACE_H
#define INTERFACE_H
#include "Structures.h"
#include <map>
#include <chrono>
#include <ctime>
#include <cstdlib>
#include <sstream>
#include <iomanip>

// inputReader.cpp
Error readParameter(std::string filename, Parameter* parameter);
Error readGeometry(std::string filename, Field* currentField, std::vector<double>* volVector);

// writeField.cpp
std::string creatDirectory(std::string dirname);

void writeField(Field* field, double t, Parameter* parameter,
                std::string const &parameterFilename="Undefined",
                std::string const &geometryFilename="Undefined",
                std::string const &filename="result");

void paraView(std::string const &filename,
              int step,
              std::vector<double> (&pos)[3],
              std::map<std::string, std::vector<double> *> const &scalars,
              std::map<std::string, std::vector<double> (*)[3] > const &vectors);

void matlab(std::string const &filename,
              std::string const &parameterFilename,
              std::string const &geometryFilename,
              int step, Parameter* parameter, Field* field);

// ConsistencyCheck.cpp
Error consistency(Parameter* param, Field* field);


#endif
