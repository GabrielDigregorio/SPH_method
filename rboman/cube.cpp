#include "sph.h"

// build a cube of particles aligned with x,y,z axes.
//  - o[3]: corner of the cube with the lowest (x,y,z) values
//  - L[3]: edge lengths along x,y and z
//  - s: particle spacing

void meshcube(double o[3], double L[3], double s, std::vector<double> &pos)
{
    // calculate nb of particles along each direction from target size "s"
    int ni = int(ceil(L[0]/s));
    double dx = L[0]/ni; ++ni;
    int nj = int(ceil(L[1]/s));
    double dy = L[1]/nj; ++nj;
    int nk = int(ceil(L[2]/s));
    double dz = L[2]/nk; ++nk;

    // output
    std::cout << "meshing cube at o=(" <<o[0]<< ","  <<o[1]<< ","  <<o[2]<< ") ";
    std::cout << "of size L=(" <<L[0]<< ","  <<L[1]<< ","  <<L[2]<< ")\n";
    std::cout << "\tparticle spacing s=(" <<dx<< ","  <<dy<< ","  <<dz<< ") [target was s=" << s << "]\n";
    std::cout << "\t=> "<<ni<< "*"  <<nj<< "*"  <<nk<< " = " << ni*nj*nk << " particles to be generated\n";

    // memory allocation
    pos.reserve(pos.size() + ni*nj*nk*3); 

    // particle generation
    for(int i=0; i<ni; ++i)
    {
        double x = o[0]+ i*dx;
        for(int j=0; j<nj; ++j)
        {
            double y = o[1]+j*dy;
            for(int k=0; k<nk; ++k)
            {
                double z = o[2]+k*dz;
                pos.push_back(x);
                pos.push_back(y);
                pos.push_back(z);
            }
        }
    }
}


