#include "sph.h"
#include <fstream>
#include <sstream>
#include <iomanip>

// export results to paraview (VTK polydata - legacy file fomat)
//   filename: file name without vtk extension
//   pos:     positions (vector of size 3*number of particles)
//   step:    time step number
//   scalars: scalar fields defined on particles (map linking [field name] <=> [vector of results v1, v2, v3, v4, ...]
//   vectors: vector fields defined on particles (map linking [field name] <=> [vector of results v1x, v1y, v1z, v2x, v2y, ...]

void paraview(std::string const &filename, 
              int step,
              std::vector<double> const &pos,
              std::map<std::string, std::vector<double> *> const &scalars,
              std::map<std::string, std::vector<double> *> const &vectors)
{
    int nbp = pos.size()/3;
    assert(pos.size()==nbp*3); // should be multiple of 3
    
    // build file name + stepno + vtk extension
    std::stringstream s; s << filename << std::setw(8) << std::setfill('0') << step << ".vtk";

    // open file
    std::cout << "writing results to " << s.str() << '\n';
    std::ofstream f(s.str().c_str());
    f << std::scientific;
    // header
    f << "# vtk DataFile Version 3.0\n";
    f << "file written by sph.exe\n";
    f << "ASCII\n";
    f << "DATASET POLYDATA\n";

    // points
    f << "POINTS " << nbp << " float\n";
    for(int i=0; i<nbp; ++i)
        f << pos[3*i+0] << " " << pos[3*i+1] << " " << pos[3*i+2] << '\n';

    // vertices
    f << "VERTICES " << nbp << " " << 2*nbp << "\n";
    for(int i=0; i<nbp; ++i)
        f << "1 " << i << '\n';
    f << '\n'; // empty line (required)

    // fields
    f << "POINT_DATA " << nbp << '\n';
    f << "FIELD FieldData " << scalars.size()+vectors.size() << '\n';

    // scalar fields
    std::map<std::string, std::vector<double> *>::const_iterator it=scalars.begin();
    for(; it!=scalars.end(); ++it)
    {
        assert(it->second->size()==nbp);
        f << it->first << " 1 " << nbp << " float\n";
        for(int i=0; i<nbp; ++i)
            f << (*it->second)[i] << '\n';
    }

    // vector fields
    it = vectors.begin();
    for(; it!=vectors.end(); ++it)
    {
        assert(it->second->size()==3*nbp);
        f << it->first << " 3 " << nbp << " float\n";
        for(int i=0; i<nbp; ++i)
            f << (*it->second)[3*i+0] << " " << (*it->second)[3*i+1] << " " << (*it->second)[3*i+2] << '\n';
    }
    f.close();
}
