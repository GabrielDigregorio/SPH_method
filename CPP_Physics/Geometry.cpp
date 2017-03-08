#include "Main.h"
#include "Physics.h"

// Build a brick of regulary aligned particles with the center of mass at o(x,y,z).
//  - o[3]: corner of the cube with the lowest (x,y,z) values ??? Center ???
//  - L[3]: edge lengths along x,y and z
//  - s: particle spacing
//  - pertubation: percentage of perturbation in position of particles (equal 0 by default)

void meshcube(double o[3], double L[3], double s, std::vector<double> &pos, double perturbation, bool stack)
{
    // open a file to write the geometry (check for valydity) MUST BE REMOVE LATER
    std::ofstream myfile;
    myfile.open ("./Results/Playground.txt", std::ofstream::out | std::ofstream::app);

    // if we stack the cube:
    if(stack == true){
        L[0] -= s/2; L[1] -= s/2; L[2] -= s/2;
    }

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
    std::default_random_engine generator; // A SEED MUST BE USE TO CHANGE VALUE AT EACH CALL 
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
                myfile << pos.end()[-3] << " " << pos.end()[-2]  << " " << pos.end()[-1]  << "\n" ;
            }
        }
    }

    myfile.close();
}



// Build a cylinder of regulary aligned particles with the center of mass at o(x,y,z).
//  - o[3]: center of the cylinder
//  - L[3]: diameter1 of the cylinder, diameter2 of the cylinder, length of the cylinder
//  - s: particle spacing
//  - (optional) pertubation: percentage of perturbation in position of particles (equal 0 by default)
//  - (optional) stack: reduce L[3] by s/2 in order to stack cube (equal 0 by default)

void meshcylinder(double o[3], double L[3], double s, std::vector<double> &pos, double perturbation, bool stack)
{
    // open a file to write the geometry (check for valydity) MUST BE REMOVE LATER
    std::ofstream myfile;
    myfile.open ("./Results/Playground.txt", std::ofstream::out | std::ofstream::app);

    // if we stack the cylinder:
    if(stack == true){
        L[0] -= s/2; L[1] -= s/2; L[2] -= s/2;
    }

    // ellipse parameter
    double a=L[0]/2, b=L[1]/2;

    // calculate nb of particles along the radius from target size "s"
    int nd1 = int(ceil(L[0]/s));
    double dr1 = L[0]/nd1; ++nd1;
    int nd2 = int(ceil(L[1]/s));
    double dr2 = L[1]/nd2; ++nd2;
    int nl = int(ceil(L[2]/s));
    double dl = L[2]/nl; ++nl;

    // output
    std::cout << "meshing cylinder at o=(" <<o[0]<< ","  <<o[1]<< ","  <<o[2]<< ") ";
    std::cout << "of diameter D1=" << L[0] << ", diameter D2=" << L[1] << " and L=" << L[2] << "\n";
    std::cout << "\tparticle spacing s=(" <<dr1<< " and "<<dr2<< "," <<dl<< ") [target was s=" << s << "]\n";
    std::cout << "\t less than  => "<<nl<< "*"  <<nd1<< "*"  <<nd2<< " = " << nl*nd1/2*nd2/2 << " particles to be generated\n";

    // memory allocation
    pos.reserve(pos.size() + nl*nd1/2*nd2/2*3);

    // generates number in the range -s*perturbation % and s*perturbation %
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(-s*perturbation/100,s*perturbation/100);

    // particle generation
    for(int l=(-nl+1)/2; l<=(nl-1)/2; ++l)
    {
        double z = o[2] + l*dl;
        for(int i=(-nd1+1)/2 ; i<=(nd1-1)/2; ++i)
        {
            double tmpx = i*dr1;
            for(int j= (-nd2+1)/2; j<=(nd2-1)/2; ++j)
            {
                if(j*s>- sqrt(pow(b,2)*(1-(pow(tmpx,2)/pow(a,2)))) && j*s<= sqrt( pow(b,2)*(1-(pow(tmpx,2)/pow(a,2)))) )
                {
                    double x = o[0] + i*dr1;
                    double y = o[1] + j*dr2;
                    pos.push_back(x + distribution(generator));
                    pos.push_back(y + distribution(generator));
                    pos.push_back(z + distribution(generator));
                    myfile << pos.end()[-3] << " " << pos.end()[-2]  << " " << pos.end()[-1]  << "\n" ;
                }
            }
        }
    }

    myfile.close();
}


// !!!DOESN't WORK!!! Build a sphere of regulary aligned particles with the center of mass at o(x,y,z).
//  - o[3]: center of the cylinder
//  - L[3]: diameter1 of the sphere, diameter2 of the sphere, , diameter3 of the sphere
//  - s: particle spacing
//  - (optional) pertubation: percentage of perturbation in position of particles (equal 0 by default)
//  - (optional) stack: reduce L[3] by s/2 in order to stack cube (equal 0 by default)

void meshsphere(double o[3], double L[3], double s, std::vector<double> &pos, double perturbation, bool stack)
{
    // open a file to write the geometry (check for valydity) MUST BE REMOVE LATER
    std::ofstream myfile;
    myfile.open ("./Results/Playground.txt", std::ofstream::out | std::ofstream::app);

    // if we stack the cylinder:
    if(stack == true){
        L[0] -= s/2; L[1] -= s/2; L[2] -= s/2;
    }

    // ellipse parameter
    double a=L[0]/2, b=L[1]/2, c=L[2]/2;

    // calculate nb of particles along the radius from target size "s"
    int nd1 = int(ceil(L[0]/s));
    double dr1 = L[0]/nd1; ++nd1;
    int nd2 = int(ceil(L[1]/s));
    double dr2 = L[1]/nd2; ++nd2;
    int nd3 = int(ceil(L[2]/s));
    double dr3 = L[2]/nd3; ++nd3;

    // output
    std::cout << "meshing sphere at o=(" <<o[0]<< ","  <<o[1]<< ","  <<o[2]<< ") ";
    std::cout << "of diameter D1=" << L[0] << ", diameter D2=" << L[1] << ", diameter D3=" << L[2] << "\n";
    std::cout << "\tparticle spacing s=(" <<dr1<< " and "<<dr2<< "," <<dr3<< ") [target was s=" << s << "]\n";
    std::cout << "\t less than => "<<nd1<< "*"  <<nd2<< "*"  <<nd3<< " = " <<nd1/2*nd2/2*nd3/2<< " particles to be generated\n";

    // memory allocation
    pos.reserve(pos.size() + nd1/2*nd2/2*nd3/2*3);

    // generates number in the range -s*perturbation % and s*perturbation %
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(-s*perturbation/100,s*perturbation/100);

    // particle generation
    for(int l=(-nd3+1)/2; l<=(nd3-1)/2; ++l)
    {
        double tmpz = l*dr3;
        double tempo=sqrt(1-pow(tmpz,2)/pow(c,2));
        double borne_y=floor((tempo*b)/(2*dr2));
        for(int j= -borne_y; j<=borne_y; ++j)
            {
                double tmpy = j*dr2;
                double tempo2=sqrt(pow(tempo,2)*pow(a,2)-pow(tmpy,2)*pow(a,2)/pow(b,2));
                double borne_x=floor(tempo2/(2*dr1));
                for(int i=-borne_x; i<=borne_x; ++i)
                {
                        double x = o[0] + i*dr1;
                        double y = o[1] + j*dr2;
                        double z = o[2] + l*dr3;
                        pos.push_back(x + distribution(generator));
                        pos.push_back(y + distribution(generator));
                        pos.push_back(z + distribution(generator));
                        myfile << pos.end()[-3] << " " << pos.end()[-2]  << " " << pos.end()[-1]  << "\n" ;
                }
            }
    }


    myfile.close();

}
