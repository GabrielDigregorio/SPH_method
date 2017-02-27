#include "SPH.hpp"


using namespace std;

// build a cube of particles aligned with x,y,z axes.
//  - o[3]: corner of the cube with the lowest (x,y,z) values ??? Center ???
//  - L[3]: edge lengths along x,y and z
//  - s: particle spacing
//  - pertubation: percentage of perturbation in position of particles (equal 0 by default)

void meshcube(double o[3], double L[3], double s, std::vector<double> &pos, double perturbation)
{
    // open a file to write the geometry (check for valydity) MUST BE REMOVE LATER
    ofstream myfile;
    myfile.open ("cube.txt");

    // calculate nb of particles along each direction from target size "s"
    int ni = int(ceil(L[0]/s));
    double dx = L[0]/ni; ++ni;
    int nj = int(ceil(L[1]/s));
    double dy = L[1]/nj; ++nj;
    int nk = int(ceil(L[2]/s));
    double dz = L[2]/nk; ++nk;

    // output
    std::cout << "meshing cube at center o=(" << o[0] << ","  << o[1] << ","  << o[2] << ") ";
    std::cout << "of size L=(" <<L[0]<< ","  <<L[1]<< ","  <<L[2]<< ")\n";
    std::cout << "\tparticle spacing s=(" <<dx<< ","  <<dy<< ","  <<dz<< ") [target was s=" << s << "]\n";
    std::cout << "\t=> "<<ni<< "*"  <<nj<< "*"  <<nk<< " = " << ni*nj*nk << " particles to be generated\n";
    
    // memory allocation
    pos.reserve(pos.size() + ni*nj*nk*3);

    // generates number in the range -s*perturbation % and s*perturbation %
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(-s*perturbation/100,s*perturbation/100);

    // particle generation
    for(int k=0; k<nk; ++k)
    {
        double z = o[2] - L[2]/2 + k*dz;
        for(int j=0; j<nj; ++j)
        {
            double y = o[1] - L[1]/2 +j*dy;
            for(int i=0; i<ni; ++i)
            {
                double x = o[0] - L[0]/2 +i*dx;
                pos.push_back(x + distribution(generator));
                pos.push_back(y + distribution(generator));
                pos.push_back(z + distribution(generator));
                //myfile << pos.end()[-3] << " " << pos.end()[-2]  << " " << pos.end()[-1]  << "\n" ;
            }
        }
    }

    myfile.close();

}



// build a cylinder of particles with the center of the lower base aligned with x,y,z axes.
//  - o[3]: center of the cylinder
//  - L[3]: diameter1 of the cylinder, diameter2 of the cylinder, length of the cylinder
//  - s: particle spacing
//  - pertubation: percentage of perturbation in position of particles (equal 0 by default)

void meshcylinder(double o[3], double L[3], double s, std::vector<double> &pos, double perturbation)
{
    // open a file to write the geometry (check for valydity) MUST BE REMOVE LATER
    ofstream myfile;
    myfile.open ("cylinder.txt");

    // calculate nb of particles along the radius from target size "s"
    int nl = int(ceil(L[2]/s));
    double dl = L[2]/nl; ++nl;
    int nr = int(ceil(L[0]/(2*s)));
    double dr = L[0]/(2*nr); ++nr;

    // output
    std::cout << "meshing cylinder at o=(" <<o[0]<< ","  <<o[1]<< ","  <<o[2]<< ") ";
    std::cout << "of diameter D1=" << L[0] << ", diameter D2=" << L[1] << " and L=" << L[2] << "\n";
    std::cout << "\tparticle spacing s=(" <<dr<< "," <<dl<< ") [target was s=" << s << "]\n";
    std::cout << "\t=> "<<nl<< "*"  <<nr<< "*"  <<nr<< " = " << nl*nr*nr << " particles to be generated\n";

    // memory allocation
    pos.reserve(pos.size() + nl*nr*nr*3);

    // generates number in the range -s*perturbation % and s*perturbation %
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(-s*perturbation/100,s*perturbation/100);

    // particle generation
    for(int l=(-nl+1)/2; l<=(nl-1)/2; ++l)
    {
        double z = o[2]+ l*dl;
        for(int i=-nr+1; i<nr; ++i)
        {
            double x = o[0]+i*dr;
            for(int j=-int(sqrt(pow(L[1]/2,2)-pow(x,2))); j<int(sqrt(pow(L[1]/2,2)-pow(x,2)))+1; ++j)
            {
                double y = o[1]+j*dr;
                pos.push_back(x + distribution(generator));
                pos.push_back(y + distribution(generator));
                pos.push_back(z + distribution(generator));
                //myfile << pos.end()[-3] << " " << pos.end()[-2]  << " " << pos.end()[-1]  << "\n" ;
            }
        }
    }

    myfile.close();

}

