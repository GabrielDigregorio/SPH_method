#include "Main.h"
#include "Interface.h"
#include "Tools.h"


// Creat directory to store data
// In: name of the directory
/* std::string creatDirectory(std::string dirname){

    std::stringstream newdir, outdir; outdir<< dirname; newdir <<"Results/"<< dirname;
    DIR* dir = opendir(newdir.str().c_str());
    int i=1;

    while(dir)
    {
         Directory exists.
        closedir(dir);
        newdir << i;
        outdir << i;
        dir = opendir(newdir.str().c_str());
        i++;
    }

    mode_t nMode = 0733; // UNIX style permissions
    int nError = 0;
    #if defined(_WIN32)
    nError = _mkdir(newdir.str().c_str()); // can be used on Windows
    #else
    nError = mkdir(newdir.str().c_str(),nMode); // can be used on non-Windows
    #endif
    // handle your error here
    }

    //mkdir(newdir.str().c_str());
    outdir<< "/"<<dirname;
    std::cout <<"\n"<<  outdir.str()<<"\n";
    return outdir.str();
} */

/*
 * In: field = stucture containing value to write
 *     t = time corresponding to the file to write
 *     filename = Name given to the file
 *     parameterFilename = Fluid parameter file used
 *     geometryFilename = geometry file used
 * Out: speed_t.vtk, pos_t.vtk, or .txt
 */
void writeField(Field* field, double t, Parameter* parameter,
                std::string const &parameterFilename,
                std::string const &geometryFilename,
                std::string const &filename)
{
    std::map<std::string, std::vector<double> *> scalars;
    std::map<std::string, std::vector<double> *> vectors;

    // Map particules
    if(parameter->format==1 || parameter->format==3)
    {
        scalars["pressure"] = &field->pressure;
        scalars["density"]  = &field->density;
        vectors["velocity"] = &field->speed;
    }

    // Save results to disk (ParaView or Matlab)
    switch (parameter->format){

    case ParaView : // .vtk in ParaView
        paraView(filename, t, field->pos, scalars, vectors);
        //return ;
    break;

    case Matlab : // .txt in Matlab
        matlab(filename, parameterFilename, geometryFilename,  t, parameter, field);
        //return ;
    break;

    case Both : // .vtk in ParaView and .txt in Matlab
        paraView(filename, t, field->pos, scalars, vectors);
        matlab(filename, parameterFilename, geometryFilename, t, parameter, field);
        //return ;
    break;

    default:
        std::cout<<"Non existing writing type."<<std::endl;
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
    std::stringstream s; s <<"Results/"<< filename << "_" << std::setw(8) << std::setfill('0') << step << ".vtk";

    // open file
    //std::cout << "writing results to " << s.str() << std::endl;
    std::ofstream f(s.str().c_str());
    f << std::scientific;
    // header
    f << "# vtk DataFile Version 3.0"<<std::endl;
    f << "file written by sph.exe"<<std::endl;
    f << "ASCII"<<std::endl;
    f << "DATASET POLYDATA"<<std::endl;

    // points
    f << "POINTS " << nbp << " float"<<std::endl;
    for(int i=0; i<nbp; ++i)
        f << pos[3*i+0] << " " << pos[3*i+1] << " " << pos[3*i+2] << std::endl;

    // vertices
    f << "VERTICES " << nbp << " " << 2*nbp << std::endl;
    for(int i=0; i<nbp; ++i)
        f << "1 " << i << std::endl;
    f << '\n'; // empty line (required)

    // fields
    f << "POINT_DATA " << nbp << '\n';
    f << "FIELD FieldData " << scalars.size()+vectors.size() << std::endl;

    // scalar fields
    std::map<std::string, std::vector<double> *>::const_iterator it=scalars.begin();
    for(; it!=scalars.end(); ++it)
    {
        assert(it->second->size()==nbp);
        f << it->first << " 1 " << nbp << " float"<<std::endl;
        for(int i=0; i<nbp; ++i)
            f << (*it->second)[i] << '\n';
    }

    // vector fields
    it = vectors.begin();
    for(; it!=vectors.end(); ++it)
    {
        assert(it->second->size()==3*nbp);
        f << it->first << " 3 " << nbp << " float"<<std::endl;
        for(int i=0; i<nbp; ++i)
            f << (*it->second)[3*i+0] << " " << (*it->second)[3*i+1] << " " << (*it->second)[3*i+2] << std::endl;
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
              int step, Parameter* parameter,Field* field)
{
    int nbp = field->pos.size()/3;
    assert(field->pos.size()==nbp*3); // should be multiple of 3

    // Set Chronos and time variable
    std::chrono::time_point<std::chrono::system_clock> start, end;

    // Start Chrono
    start = std::chrono::system_clock::now();

    // Date
    std::time_t Date = std::chrono::system_clock::to_time_t(start);
    time_t rawtime; struct tm * timeinfo; time (&rawtime); timeinfo = localtime (&rawtime);

    // build file name + stepno + vtk extension
    std::stringstream s; s << "Results/"<< filename << "_" << std::setw(8) << std::setfill('0') << step << ".txt";

    // open file
    //std::cout << "Writing results to " << s.str() << std::endl;
    std::ofstream f(s.str().c_str());
    f << std::scientific;

    //Record Time
    double duration = ( std::clock() - startExperimentTimeClock ) / (double) CLOCKS_PER_SEC;
    // header
    f << "#EXPERIMENT: " << filename << std::endl;
    f << std::endl;
    f << "Date : " << asctime(timeinfo);
    #if defined(_WIN32) || defined(WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
    f << "Computer Name : "<< getenv("COMPUTERNAME") <<std::endl;
    f << "Username : "<< getenv("USERNAME") <<std::endl;
    #elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
    f << "Computer Name : None"<<std::endl;// To be implemented
    f << "Username : None"<<std::endl;// To be implemented
    #else
    #error "Cannot define GetMemory( ) or GetMemoryProcessPeak( ) or GetMemoryProcess() for an unknown OS."
    #endif
    f << "File Used : " << geometryFilename << "   &   " << parameterFilename << std::endl;
    f << std::endl;
    f << "CPU Time : " << duration << " [s]" << std::endl;
    f << "Memory Usage : " << GetMemoryProcess(false, false)<< " [kB]" <<std::endl;
    f << "Memory Usage Peak : " << GetMemoryProcessPeak(false, false)<< " [kB]" <<std::endl;
    f << std::endl;
    f << "Step Time (k) : "<< parameter->k << " [s]" << std::endl;
    f << "Write interval : "<< parameter->writeInterval << " [s]" << std::endl;
    f << "Simulation Time (T) : "<< parameter->T << " [s]" << std::endl;
    f << std::endl;
    f << "Domain (lower l) : "<< field->l[0] << "   " << field->l[1] << "   " << field->l[2] << "    [m]" << std::endl;
    f << "Domain (upper u) : "<< field->u[0] << "   " << field->u[1] << "   " << field->u[2] << "    [m]" << std::endl;
    f << "Number of Particules (nFree/nMoving/nFixed) : "<< field->nFree << "   "<< field->nMoving<< "   "<< field->nFixed <<  std::endl;
    f << "\n";
    f << " posX\t        posY\t        posZ\t     velocityX\t     velocityY\t     velocityZ\t     density\t     pressure\t     mass"<<std::endl;

    // Fill f:
    for(int i=0; i<nbp; ++i)
    {
        f <<field->pos[3*i+0]<<"\t"<<field->pos[3*i+1]<<"\t"<<field->pos[3*i+2]<<"\t"
        <<field->speed[3*i+0]<<"\t"<<field->speed[3*i+1]<<"\t"<<field->speed[3*i+2]<<"\t"
        <<field->density[i]<<"\t"
        <<field->pressure[i]<<"\t"
        <<field->mass[i]<<"\t"<<std::endl;
    }


    // End Chrono
    end = std::chrono::system_clock::now();

    // Write result on file
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);
    //std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";

    f.close();

    }
