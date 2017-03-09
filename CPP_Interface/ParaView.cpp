#include "Main.h"
#include "Interface.h"


// export results to paraview (VTK polydata - legacy file fomat)
//   filename: file name without vtk extension
//   pos:     positions (vector of size 3*number of particles)
//   t:       time number
//   scalars: scalar fields defined on particles (map linking [field name] <=> [vector of results v1, v2, v3, v4, ...]
//   vectors: vector fields defined on particles (map linking [field name] <=> [vector of results v1x, v1y, v1z, v2x, v2y, ...]

void paraview(std::string const &filename, 
              double t,
              std::vector<double> const &posFree,
              std::map<std::string, std::vector<double> *> const &scalarsFree,
              std::map<std::string, std::vector<double> *> const &vectorsFree)
{
    int nbpFree = posFree.size()/3;
    assert(posFree.size()==nbpFree*3); // should be multiple of 3
    
    // build file name + time_no + vtk extension
    std::stringstream s; s << filename << std::setw(8) << std::setfill('0') << t << ".vtk";

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
    f << "POINTS " << nbpFree << " float\n";
    for(int i=0; i<nbpFree; ++i)
        f << posFree[3*i+0] << " " << posFree[3*i+1] << " " << posFree[3*i+2] << '\n';

    // vertices
    f << "VERTICES " << nbpFree << " " << 2*nbpFree << "\n";
    for(int i=0; i<nbp; ++i)
        f << "1 " << i << '\n';
    f << '\n'; // empty line (required)

    // fields
    f << "POINT_DATA " << nbpFree << '\n';
    f << "FIELD FieldData " << scalarsFree.size()+vectorsFree.size() << '\n';

    // scalar fields
    std::map<std::string, std::vector<double> *>::const_iterator it=scalarsFree.begin();
    for(; it!=scalarsFree.end(); ++it)
    {
        assert(it->second->size()==nbpFree);
        f << it->first << " 1 " << nbpFree << " float\n";
        for(int i=0; i<nbpFree; ++i)
            f << (*it->second)[i] << '\n';
    }

    // vector fields
    it = vectorsFree.begin();
    for(; it!=vectorsFree.end(); ++it)
    {
        assert(it->second->size()==3*nbpFree);
        f << it->first << " 3 " << nbpFree << " float\n";
        for(int i=0; i<nbpFree; ++i)
            f << (*it->second)[3*i+0] << " " << (*it->second)[3*i+1] << " " << (*it->second)[3*i+2] << '\n';
    }
    f.close();
}
