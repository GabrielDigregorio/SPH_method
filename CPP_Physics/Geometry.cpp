#include "Main.h"
#include "Physics.h"
#define M_PI           3.14159265358979323846  /* pi */
// Build a brick of regulary aligned particles with the center of mass at o(x,y,z).
//  - o[3]: corner of the cube with the lowest (x,y,z) values ??? Center ???
//  - L[3]: edge lengths along x,y and z
//  - s: particle spacing
//  - pertubation: percentage of perturbation in position of particles (equal 0 by default)

/*void meshcube(double o[3], double L[3], double s, std::vector<double> &pos, int* nPart, double* volPart,
     double perturbation, bool stack){
    // if we stack the cube:
    if(stack == true){
        L[0] -= s; L[1] -= s; L[2] -= s;
    }

    // calculate nb of particles along each direction from target size "s"
    int ni = int(ceil(L[0]/s));
    double dx = L[0]/ni; ++ni;
    int nj = int(ceil(L[1]/s));
    double dy = L[1]/nj; ++nj;
    int nk = int(ceil(L[2]/s));
    double dz = L[2]/nk; ++nk;
    // Volume & number of particles computation
    (*nPart)=ni*nj*nk;
    (*volPart)=dx*dy*dz;


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
            }
        }
    }
}*/



// Build a cylinder of regulary aligned particles with the center of mass at o(x,y,z).
//  - o[3]: center of the cylinder
//  - L[3]: diameter1 of the cylinder, diameter2 of the cylinder, length of the cylinder
//  - s: particle spacing
//  - (optional) pertubation: percentage of perturbation in position of particles (equal 0 by default)
//  - (optional) stack: reduce L[3] by s/2 in order to stack cube (equal 0 by default)

void meshcylinder(double o[3], double L[3], double s, std::vector<double> &pos, int* nPart, double* volPart,
     double perturbation, bool stack){
    // if we stack the cylinder:
    if(stack == true){
        L[0] -= s; L[1] -= s; L[2] -= s;
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
    // Volume & number of particles computation
    (*nPart)=nl*nd1/2*nd2/2;
    //(*volPart)=;

    /* output
    std::cout << "meshing cylinder at o=(" <<o[0]<< ","  <<o[1]<< ","  <<o[2]<< ") ";
    std::cout << "of diameter D1=" << L[0] << ", diameter D2=" << L[1] << " and L=" << L[2] << "\n";
    std::cout << "\tparticle spacing s=(" <<dr1<< " and "<<dr2<< "," <<dl<< ") [target was s=" << s << "]\n";
    std::cout << "\t less than  => "<<nl<< "*"  <<nd1<< "*"  <<nd2<< " = " << (*nPart) << " particles to be generated\n";
    */

    // memory allocation
    pos.reserve(pos.size() + (*nPart)*3);

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
                }
            }
        }
    }
}


// !!!DOESN't WORK!!! Build a sphere of regulary aligned particles with the center of mass at o(x,y,z).
//  - o[3]: center of the cylinder
//  - L[3]: diameter1 of the sphere, diameter2 of the sphere, , diameter3 of the sphere
//  - s: particle spacing
//  - (optional) pertubation: percentage of perturbation in position of particles (equal 0 by default)
//  - (optional) stack: reduce L[3] by s/2 in order to stack cube (equal 0 by default)

void meshsphere(double o[3], double L[3], double s, std::vector<double> &pos, int* nPart, double* volPart,
     double perturbation, bool stack){
           // if we stack the cylinder:
    if(stack == true){
        L[0] -= s; L[1] -= s; L[2] -= s;
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
    // Volume & number of particles computation
    (*nPart)=nd1/2*nd2/2*nd3/2;
    //(*volPart)=;

    /* output
    std::cout << "meshing sphere at o=(" <<o[0]<< ","  <<o[1]<< ","  <<o[2]<< ") ";
    std::cout << "of diameter D1=" << L[0] << ", diameter D2=" << L[1] << ", diameter D3=" << L[2] << "\n";
    std::cout << "\tparticle spacing s=(" <<dr1<< " and "<<dr2<< "," <<dr3<< ") [target was s=" << s << "]\n";
    std::cout << "\t less than => "<<nd1<< "*"  <<nd2<< "*"  <<nd3<< " = " <<(*nPart)<< " particles to be generated\n";
    */

    // memory allocation
    pos.reserve(pos.size() + (*nPart)*3);

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
                }
            }
    }
}
void meshcube(double o[3], double L[3],double teta[3],double s, std::vector<double> &pos, int* nPart, double* volPart, double perturbation, bool stack){
     // the function is the same as meshcube, but as an option to make the palne rotate in any direction over it's center à masse
      // if we stack the cube:

      //std::cout<<"teta " << teta[0]<< teta[1]<< teta[2] <<"\n";
    if(stack == true){
        L[0] -= s; L[1] -= s; L[2] -= s;
    }
    int flag_1 =0;
    int flag_2 =0;
    int flag_3 =0;
   if(L[0]==0)
   {
   L[0]=s;
   flag_1=1;
   }
   if(L[1]==0){
   L[1]=s;
   flag_2=1;
   }
   if(L[2]==0)
   {
   L[2]=s;
   flag_3=1;
   }


    // calculate nb of particles along each direction from target size "s"

    int ni = int(ceil(L[0]/s));
    double dx = L[0]/ni; ++ni;
    if(flag_1==1){
    ni=ni-1;}
    int nj = int(ceil(L[1]/s));
    double dy = L[1]/nj; ++nj;
    if(flag_2==1){
    nj=nj-1;}
    int nk = int(ceil(L[2]/s));
    double dz = L[2]/nk; ++nk;
    if(flag_3==1){
    nk=nk-1;}
    // Volume & number of particles computation
    (*nPart)=ni*nj*nk;
    (*volPart)=dx*dy*dz;

    // output
    /*
    std::cout << "meshing cube at center o=(" << o[0] << ","  << o[1] << ","  << o[2] << ") ";
    std::cout << "of size L=(" <<L[0]<< ","  <<L[1]<< ","  <<L[2]<< ")\n";
    std::cout << "\tparticle spacing s=(" <<dx<< ","  <<dy<< ","  <<dz<< ") [target was s=" << s << "]\n";
    std::cout << "\t=> "<<ni<< "*"  <<nj<< "*"  <<nk<< " = " << (*nPart) << " particles to be generated\n";
    */

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
            }
        }
    }
    int n_block=3*nk*nj*ni;
    int n_total=pos.size();
    int start=(n_total-n_block)/3;
    // creation matrice rotation
    std::vector<double> Rx(9,0);
    std::vector<double> Ry(9,0);
    std::vector<double> Rz(9,0);
    double teta_x=teta[0];// en dregré
    double teta_y=teta[1];// en dregré
    double teta_z=teta[2];// en dregré
    Rx[0]=1;
    Rx[1]=0;
    Rx[2]=0;
    Rx[3]=0;
    Rx[4]=cos(teta_x/180.0*M_PI);
    Rx[5]=-sin(teta_x/180.0*M_PI);
    Rx[6]=0;
    Rx[7]=sin(teta_x/180.0*M_PI);
    Rx[8]=cos(teta_x/180.0*M_PI);
    Ry[0]=cos(teta_y/180.0*M_PI);
    Ry[1]=0;
    Ry[2]=sin(teta_y/180.0*M_PI);
    Ry[3]=0;
    Ry[4]=1;
    Ry[5]=0;
    Ry[6]=-sin(teta_y/180.0*M_PI);
    Ry[7]=0;
    Ry[8]=cos(teta_y/180.0*M_PI);

    Rz[0]=cos(teta_z/180.0*M_PI);
    Rz[1]=-sin(teta_z/180.0*M_PI);
    Rz[2]=0;
    Rz[3]=sin(teta_z/180.0*M_PI);
    Rz[4]=cos(teta_z/180.0*M_PI);
    Rz[5]=0;
    Rz[6]=0;
    Rz[7]=0;
    Rz[8]=1;
        for(int i=start; i<(n_total)/3; ++i)
            {
                pos[3*i]= pos[3*i]-o[0];
                pos[3*i+1]= pos[3*i+1]-o[1];
                pos[3*i+2]= pos[3*i+2]-o[2];
            }
        // rotation
        for(int i=start; i<(n_total)/3; ++i)
        {
           // first Rx
            double var_x=pos[3*i];
            double var_y=pos[3*i+1];
            double var_z=pos[3*i+2];
            pos[3*i]= Rx[0]*var_x+Rx[1]*var_y+Rx[2]*var_z;
            pos[3*i+1]= Rx[3]*var_x+Rx[4]*var_y+Rx[5]*var_z;
            pos[3*i+2]= Rx[6]*var_x+Rx[7]*var_y+Rx[8]*pos[3*i+2];
            
            // second Ry
            var_x=pos[3*i];
            var_y=pos[3*i+1];
            var_z=pos[3*i+2];
            pos[3*i]= Ry[0]*var_x+Ry[1]*var_y+Ry[2]*var_z;
            pos[3*i+1]= Ry[3]*var_x+Ry[4]*var_y+Ry[5]*var_z;
            pos[3*i+2]= Ry[6]*var_x+Ry[7]*var_y+Ry[8]*var_z;
            // thrid Rz
            var_x=pos[3*i];
            var_y=pos[3*i+1];
            var_z=pos[3*i+2];
            pos[3*i]= Rz[0]*var_x+Rz[1]*var_y+Rz[2]*var_z;
            pos[3*i+1]= Rz[3]*var_x+Rz[4]*var_y+Rz[5]*var_z;
            pos[3*i+2]= Rz[6]*var_x+Rz[7]*var_y+Rz[8]*var_z;
        }
        // translation back at the first place
        for(int i=start; i<(n_total)/3; ++i)
            {
                pos[3*i]= pos[3*i]+o[0];
                pos[3*i+1]= pos[3*i+1]+o[1];
                pos[3*i+2]= pos[3*i+2]+o[2];
            }


}
