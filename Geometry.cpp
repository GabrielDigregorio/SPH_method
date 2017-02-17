#include "sph.h"
#include <fstream>
using namespace std;

// build a cube of particles aligned with x,y,z axes.
//  - o[3]: corner of the cube with the lowest (x,y,z) values
//  - L[3]: edge lengths along x,y and z
//  - s: particle spacing

void meshcube(double o[3], double L[3], double s, std::vector<double> &pos)
{   
    // open a file to write the geometry (check for valydity)
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
                myfile << x << " " << y << " " << z << ";" ; 
            }
        }
    }

}



// build a cylinder of particles with the center of the lower base aligned with x,y,z axes.
//  - o[3]: center of the cylinder
//  - R: radius of the cylinder
//  - L: length of the cylinder
//  - s: particle spacing

void meshcylinder(double o[3], double L, double R, double s, std::vector<double> &pos)
{
    // open a file to write the geometry (check for valydity)
    ofstream myfile;
    myfile.open ("cylinder.txt");

    // calculate nb of particles along the radius from target size "s"
    int nl = int(ceil(L/s));
    double dl = L/nl; ++nl;
    int nr = int(ceil(R/s));
    double dr = R/nr; ++nr;

    // output
    std::cout << "meshing cylinder at o=(" <<o[0]<< ","  <<o[1]<< ","  <<o[2]<< ") ";
    std::cout << "of radius R=" <<R<< " and L=" <<L<< "\n";
    std::cout << "\tparticle spacing s=(" <<dr<< "," <<dl<< ") [target was s=" << s << "]\n";
    std::cout << "\t=> "<<nl<< "*"  <<nr<< "*"  <<nr<< " = " << nl*nr*nr << " particles to be generated\n";

    // memory allocation
    pos.reserve(pos.size() + nl*nr*nr*3);
  
    // particle generation
    for(int l=0; l<nl-1; ++l)
    {
        double z = o[2]+ l*dl;
        for(int i=-nr+1; i<nr; ++i)
        {
            double x = o[0]+i*dr;
            for(int j=-int(sqrt(pow(R,2)-pow(x,2))); j<int(sqrt(pow(R,2)-pow(x,2)))+1; ++j)
            {
                double y = o[1]+j*dr;
                pos.push_back(x);
                pos.push_back(y);
                pos.push_back(z);
                myfile << x << " " << y << " " << z << "\n" ; 
            }
        }
    }

    myfile.close();
    
}



// build a sphere of particles centered in [x,y,z]=(0,0,0).
//  - o[3]: center of the sphere
//  - r: radius of the sphere
//  - s: particle spacing

void meshsphere(double o[3], double r, double s, std::vector<double> &pos)
{
    // calculate nb of particles along the radius from target size "s"

    // output

    // memory allocation

    // particle generation
    
}