#ifndef PARAVIEW_H
#define PARAVIEW_H

#include <string>
#include <vector>
#include <map>

enum PFormat
{
    LEGACY_TXT = 0,
    LEGACY_BIN = 1, 
    XML_BIN    = 2,
    XML_BINZ   = 3
};

void paraview(std::string const &filename, 
              int step,
              std::vector<double> const &pos,
              std::map<std::string, std::vector<double> *> const &scalars,
              std::map<std::string, std::vector<double> *> const &vectors,
              PFormat format);       

#endif // PARAVIEW_H
