#include "sph.h"

int main(int argc, char *argv[])
{
    // creation of a cube of particles
    std::vector<double> pos;
    double o[3] = { 10.0, 10.0, 10.0 };
    double L[3] = { 1.0, 2.0, 3.0 };
    double s = 0.2;

    meshcube(o, L, s, pos);

    // creation of dummy pressure/density/velocity fields &
    int nbp = pos.size()/3;
    std::vector<double> pressure(nbp);
    std::vector<double> density(nbp);
    std::vector<double> velocity(nbp*3);

    std::map<std::string, std::vector<double> *> scalars;
    std::map<std::string, std::vector<double> *> vectors;
    scalars["pressure"] = &pressure;
    scalars["density"]  = &density;
    vectors["velocity"] = &velocity;

    // time step loop
    for(int nstep=0; nstep<20; ++nstep)
    {
        double a = nstep/20.0*8*atan(1.0);

        // generate dummy results
        for(int i=0; i<nbp; ++i)
        {
            pos[3*i+2] -= 0.5/20;

            double x = pos[3*i+0];
            double y = pos[3*i+1];
            double z = pos[3*i+2];

            pressure[i] = z/2+sin(a);
            density[i] = y+cos(a)+2;
            velocity[3*i+0] = x-(o[0]+L[0]/2+L[0]*sin(a));
            velocity[3*i+1] = y-(o[1]+L[1]/2+L[1]*cos(a));
            velocity[3*i+2] = z-(o[2]+L[2]/2);
        }

        // save results to disk
        paraview("results", nstep, pos, scalars, vectors);
    }

    return 0;
}
