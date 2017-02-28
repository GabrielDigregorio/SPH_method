#include "SPH.hpp"

int main(int argc, char *argv[])
{
    // creation of a cube of particles
    bool stack = true;
    std::vector<double> pos;
    double o[3] = { 0.0, 0.0, 0.0 };
    double L[3] = { 2.0, 3.0, 4.0 };
    double s = 1;
    meshcube(o, L, s, pos, 0.0, stack);

    /*
    double o[3] = { 0.0, 0.0, 0.0 };
    double L = 10;
	  double R = 10;
    double s = 1;
    meshcylinder(o, L, R, s, pos);
    */

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
    for(int nstep=0; nstep<1; ++nstep)
    {
        double a = nstep/20.0*8*atan(1.0);

        /*
        // generate dummy results
        for(int i=0; i<nbp; ++i)
        {
            pos[3*i+2] -= 0.5/20;

            double x = pos[3*i+0];
            double y = pos[3*i+1];
            double z = pos[3*i+2];

            pressure[i] = z/2+sin(a);
            density[i] = y+cos(a)+2;
            velocity[3*i+0] = 0;
            velocity[3*i+1] = 0;
            velocity[3*i+2] = 0;
        }
        */
        double l[3] = {o[0], o[1], o[2]};
        double u[3] = {o[0]+L[0], o[1]+L[1], o[2]+L[2]};
        double kh = 1.5;
        std::vector<double> values;
        std::vector<int> row;
        std::vector<int> column;
        neighborLinkedList(pos, l, u, kh, values, row, column);

        // save results to disk
        paraview("results", nstep, pos, scalars, vectors);
    }

    return 0;
}
