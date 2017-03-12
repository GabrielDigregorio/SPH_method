#include "Main.h"
#include "Interface.h"
#include "Tools.h"

/*
 * In: field = stucture containing value to write
 *     t = time corresponding to the file to write
 *     filename = Name given to the file
 *     parameterFilename = Fluid parameter file used
 *     geometryFilename = geometry file used
 * Out: speed_t.vtk, pos_t.vtk, or .txt
 */
void writeField(Field* field, double t, Format myFormat, 
                std::string const &parameterFilename,
                std::string const &geometryFilename,
                std::string const &filename)
{

    std::map<std::string, std::vector<double> *> scalars;
    std::map<std::string, std::vector<double> *> vectors;

    // Map particules
    if(myFormat==1 || myFormat==3)
    {
        scalars["pressure"] = &field->pressure;
        scalars["density"]  = &field->density;
        vectors["velocity"] = &field->speed;
    }

    // Save results to disk (ParaView or Matlab)
    switch (myFormat){

    case 1 : // .vtk in ParaView
        paraView(filename, t, field->pos, scalars, vectors);
        //return ;
    break;

    case 2 : // .txt in Matlab
        matlab(filename, parameterFilename, geometryFilename,  t,
               field->pos, field->speed, field->density, field->pressure, field->mass);
        //return ;
    break;

    case 3 : // .vtk in ParaView and .txt in Matlab
        paraView(filename, t, field->pos, scalars, vectors);
        matlab(filename, parameterFilename, geometryFilename, t,
               field->pos, field->speed, field->density, field->pressure, field->mass);
        //return ;
    break;

    default: ;
    }
}



// export results to paraview (VTK polydata - legacy file fomat)
//   filename: file name without vtk extension
//   pos:     positions (vector of size 3*number of particles)
//   step:    time step number
//   scalars: scalar fields defined on particles (map linking [field name] <=> [vector of results v1, v2, v3, v4, ...]
//   vectors: vector fields defined on particles (map linking [field name] <=> [vector of results v1x, v1y, v1z, v2x, v2y, ...]
void paraView(std::string const &filename, 
              int step,
              std::vector<double> const &pos,
              std::map<std::string, std::vector<double> *> const &scalars,
              std::map<std::string, std::vector<double> *> const &vectors)
{
    int nbp = pos.size()/3;
    assert(pos.size()==nbp*3); // should be multiple of 3
    
    // build file name + stepno + vtk extension
    std::stringstream s; s <<  "../Results/" << filename << "_" << std::setw(8) << std::setfill('0') << step << ".vtk";

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



// export results to Matlab (.txt)
//   filename: file name without txt extension
//   pos:     positions (vector of size 3*number of particles)
//   speed:   velocity  (vector of size 3*number of particles)
//   density: density   (vector of size number of particles)
//   pressure:pressure  (vector of size number of particles)
//   mass:    mass      (vector of size number of particles)
//   step:    time step number
void matlab(std::string const &filename,
              std::string const &parameterFilename,
              std::string const &geometryFilename,
              int step,
              std::vector<double> const &pos,
              std::vector<double> const &speed,
              std::vector<double> const &density,
              std::vector<double> const &pressure,
              std::vector<double> const &mass)
{
    int nbp = pos.size()/3;
    assert(pos.size()==nbp*3); // should be multiple of 3

    // Set Chronos and time variable
    std::chrono::time_point<std::chrono::system_clock> start, end;

    // Start Chrono
    start = std::chrono::system_clock::now();

    // Date
    std::time_t Date = std::chrono::system_clock::to_time_t(start);
    time_t rawtime; struct tm * timeinfo; time (&rawtime); timeinfo = localtime (&rawtime);

    // build file name + stepno + vtk extension
    std::stringstream s; s << "../Results/" << filename << "_" << std::setw(8) << std::setfill('0') << step << ".txt";

    // open file
    std::cout << "Writing results to " << s.str() << '\n';
    std::ofstream f(s.str().c_str());
    f << std::scientific;

    // header
    f << "#EXPERIMENT: " << filename << "\n";
    f << "Date : " << asctime(timeinfo);
    #if defined(_WIN32) || defined(WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
    f << "Computer Name : "<< getenv("COMPUTERNAME") <<"\n";
    f << "Username : "<< getenv("USERNAME") <<"\n";
    #elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
    f << "Computer Name : "<< getenv("HOSTNAME") <<"\n";
    f << "Username : "<< getenv("USER") <<"\n";
    #else
    #error "Cannot define GetMemory( ) or GetMemoryProcessPeak( ) or GetMemoryProcess() for an unknown OS."
    #endif
    f << "File Used : " << geometryFilename << " & " << parameterFilename << "\n";
    f << "CPU Time : " << "- [s]" << "\n";
    f << "Memory Usage : " << GetMemoryProcess(false, false)<< " [kB]" <<"\n";
    f << "Memory Usage Peak : " << GetMemoryProcessPeak(false, false)<< " [kB]" <<"\n";
    f << "\n";
    f << " posX\t\t\t posY\t\t\t posZ\t\t\t velocityX\t\t velocityY\t\t velocityZ\t\t density\t\t pressure\t\t mass\n";

    // Fill f:
    for(int i=0; i<nbp; ++i)
    {
        f <<pos[3*i+0]<<"\t"<<pos[3*i+1]<<"\t"<<pos[3*i+2]<<"\t"
        <<speed[3*i+0]<<"\t"<<speed[3*i+1]<<"\t"<<speed[3*i+2]<<"\t"
        <<density[i]<<"\t"
        <<pressure[i]<<"\t"
        <<mass[i]<<"\t"<<"\n";
    }


    // End Chrono
    end = std::chrono::system_clock::now();
 
    // Write result on file
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);
    //std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";

    f.close();
    }
