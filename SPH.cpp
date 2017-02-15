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

    // save results to disk
    /*paraview("results", nstep, pos, scalars, vectors);*/

    return 0;
}
