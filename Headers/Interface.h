///**************************************************************************
/// HEADER: Function That Load/Read Input Files And Generate Output Files
///**************************************************************************

#ifndef INTERFACE_H
#define INTERFACE_H
#include "Structures.h"


// Playground.cpp
void readParameter(char* parameterFilename, Parameter* parameter);
void readGeometry(char* geometryFilename, Field* field);

// writeField.cpp
void writeField(Field* field,double t);


// ParaView.cpp
#include <map>
void paraview(std::string const &filename,
              int step,
              std::vector<double> const &pos,
              std::map<std::string, std::vector<double> *> const &scalars,
              std::map<std::string, std::vector<double> *> const &vectors);

#endif