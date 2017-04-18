#include "Main.h"
#include "Physics.h"
#define M_PI     3.14159265358979323846

void RotateVector(std::vector<double> &pos,double teta[3],int i )
{
    std::vector<double> Rx(9,0);
    std::vector<double> Ry(9,0);
    std::vector<double> Rz(9,0);
    double tetax=teta[0];// en dregré
    double tetay=teta[1];// en dregré
    double tetaz=teta[2];// en dregré
    Rx[0]=1;
    Rx[1]=0;
    Rx[2]=0;
    Rx[3]=0;
    Rx[4]=cos(tetax/180.0*M_PI);
    Rx[5]=-sin(tetax/180.0*M_PI);
    Rx[6]=0;
    Rx[7]=sin(tetax/180.0*M_PI);
    Rx[8]=cos(tetax/180.0*M_PI);

    Ry[0]=cos(tetay/180.0*M_PI);
    Ry[1]=0;
    Ry[2]=sin(tetay/180.0*M_PI);
    Ry[3]=0;
    Ry[4]=1;
    Ry[5]=0;
    Ry[6]=-sin(tetay/180.0*M_PI);
    Ry[7]=0;
    Ry[8]=cos(tetay/180.0*M_PI);

    Rz[0]=cos(tetaz/180.0*M_PI);
    Rz[1]=-sin(tetaz/180.0*M_PI);
    Rz[2]=0;
    Rz[3]=sin(tetaz/180.0*M_PI);
    Rz[4]=cos(tetaz/180.0*M_PI);
    Rz[5]=0;
    Rz[6]=0;
    Rz[7]=0;
    double R11=(Rz[0]*Ry[0]+Rz[1]*Ry[3]+Rz[2]*Ry[6])*Rx[0]+(Rz[0]*Ry[1]+Rz[1]*Ry[4]+Rz[2]*Ry[7])*Rx[3]+(Rz[0]*Ry[2]+Rz[1]*Ry[5]+Rz[2]*Ry[8])*Rx[6];
    double R21=(Rz[3]*Ry[0]+Rz[4]*Ry[3]+Rz[5]*Ry[6])*Rx[0]+(Rz[3]*Ry[1]+Rz[4]*Ry[4]+Rz[5]*Ry[7])*Rx[3]+(Rz[3]*Ry[2]+Rz[4]*Ry[5]+Rz[5]*Ry[8])*Rx[6];
    double R31=(Rz[6]*Ry[0]+Rz[7]*Ry[3]+Rz[8]*Ry[6])*Rx[0]+(Rz[6]*Ry[1]+Rz[7]*Ry[4]+Rz[8]*Ry[7])*Rx[3]+(Rz[6]*Ry[2]+Rz[7]*Ry[5]+Rz[8]*Ry[8])*Rx[6];

    double R12=(Rz[0]*Ry[0]+Rz[1]*Ry[3]+Rz[2]*Ry[6])*Rx[1]+(Rz[0]*Ry[1]+Rz[1]*Ry[4]+Rz[2]*Ry[7])*Rx[4]+(Rz[0]*Ry[2]+Rz[1]*Ry[5]+Rz[2]*Ry[8])*Rx[7];
    double R22=(Rz[3]*Ry[0]+Rz[4]*Ry[3]+Rz[5]*Ry[6])*Rx[1]+(Rz[3]*Ry[1]+Rz[4]*Ry[4]+Rz[5]*Ry[7])*Rx[4]+(Rz[3]*Ry[2]+Rz[4]*Ry[5]+Rz[5]*Ry[8])*Rx[7];
    double R32=(Rz[6]*Ry[0]+Rz[7]*Ry[3]+Rz[8]*Ry[6])*Rx[1]+(Rz[6]*Ry[1]+Rz[7]*Ry[4]+Rz[8]*Ry[7])*Rx[4]+(Rz[6]*Ry[2]+Rz[7]*Ry[5]+Rz[8]*Ry[8])*Rx[7];

    double R13=(Rz[0]*Ry[0]+Rz[1]*Ry[3]+Rz[2]*Ry[6])*Rx[2]+(Rz[0]*Ry[1]+Rz[1]*Ry[4]+Rz[2]*Ry[7])*Rx[5]+(Rz[0]*Ry[2]+Rz[1]*Ry[5]+Rz[2]*Ry[8])*Rx[8];
    double R23=(Rz[3]*Ry[0]+Rz[4]*Ry[3]+Rz[5]*Ry[6])*Rx[2]+(Rz[3]*Ry[1]+Rz[4]*Ry[4]+Rz[5]*Ry[7])*Rx[5]+(Rz[3]*Ry[2]+Rz[4]*Ry[5]+Rz[5]*Ry[8])*Rx[8];
    double R33=(Rz[6]*Ry[0]+Rz[7]*Ry[3]+Rz[8]*Ry[6])*Rx[2]+(Rz[6]*Ry[1]+Rz[7]*Ry[4]+Rz[8]*Ry[7])*Rx[5]+(Rz[6]*Ry[2]+Rz[7]*Ry[5]+Rz[8]*Ry[8])*Rx[8];

    double temp1 = R11*pos[3*i]+R12*pos[3*i+1]+R13*pos[3*i+2];
    double temp2 = R21*pos[3*i]+R22*pos[3*i+1]+R23*pos[3*i+2];
    double temp3 = R31*pos[3*i]+R32*pos[3*i+1]+R33*pos[3*i+2];
    pos[3*i] = temp1;
    pos[3*i+1] = temp2;
    pos[3*i+2] = temp3;
}

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

void meshcube(double o[3], double L[3],double teta[3],double s, std::vector<double> &pos, int* nPart, double* volPart, double perturbation, bool stack)
{
 // if we stack the cube:
 if(stack == true)
 {
     L[0] -= s; L[1] -= s; L[2] -= s;
 }
 int flag1=0;
 int flag2=0;
 int flag3=0;
 // calculate nb of particles along each direction from target size "s"
 if(L[0]==0)
 {
     L[0]=s;
     flag1=1;
 }
 int ni = int(ceil(L[0]/s));
 double dx = L[0]/ni; ++ni;
 if(flag1==1){ni=1;}
 
 if(L[1]==0)
 {
     L[1]=s;
     flag2=1;
 }
 int nj = int(ceil(L[1]/s));
 double dy = L[1]/nj; ++nj;
 if(flag2==1){nj=1;}
 
 if(L[2]==0)
 {
     L[2]=s;
     flag3=1;
 }
 int nk = int(ceil(L[2]/s));
 double dz = L[2]/nk; ++nk;
 if(flag3==1){nk=1;}
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

    // Could be much eaiser by precomputing matrix product !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
    double R11=(Rz[0]*Ry[0]+Rz[1]*Ry[3]+Rz[2]*Ry[6])*Rx[0]+(Rz[0]*Ry[1]+Rz[1]*Ry[4]+Rz[2]*Ry[7])*Rx[3]+(Rz[0]*Ry[2]+Rz[1]*Ry[5]+Rz[2]*Ry[8])*Rx[6];
    double R21=(Rz[3]*Ry[0]+Rz[4]*Ry[3]+Rz[5]*Ry[6])*Rx[0]+(Rz[3]*Ry[1]+Rz[4]*Ry[4]+Rz[5]*Ry[7])*Rx[3]+(Rz[3]*Ry[2]+Rz[4]*Ry[5]+Rz[5]*Ry[8])*Rx[6];
    double R31=(Rz[6]*Ry[0]+Rz[7]*Ry[3]+Rz[8]*Ry[6])*Rx[0]+(Rz[6]*Ry[1]+Rz[7]*Ry[4]+Rz[8]*Ry[7])*Rx[3]+(Rz[6]*Ry[2]+Rz[7]*Ry[5]+Rz[8]*Ry[8])*Rx[6];

    double R12=(Rz[0]*Ry[0]+Rz[1]*Ry[3]+Rz[2]*Ry[6])*Rx[1]+(Rz[0]*Ry[1]+Rz[1]*Ry[4]+Rz[2]*Ry[7])*Rx[4]+(Rz[0]*Ry[2]+Rz[1]*Ry[5]+Rz[2]*Ry[8])*Rx[7];
    double R22=(Rz[3]*Ry[0]+Rz[4]*Ry[3]+Rz[5]*Ry[6])*Rx[1]+(Rz[3]*Ry[1]+Rz[4]*Ry[4]+Rz[5]*Ry[7])*Rx[4]+(Rz[3]*Ry[2]+Rz[4]*Ry[5]+Rz[5]*Ry[8])*Rx[7];
    double R32=(Rz[6]*Ry[0]+Rz[7]*Ry[3]+Rz[8]*Ry[6])*Rx[1]+(Rz[6]*Ry[1]+Rz[7]*Ry[4]+Rz[8]*Ry[7])*Rx[4]+(Rz[6]*Ry[2]+Rz[7]*Ry[5]+Rz[8]*Ry[8])*Rx[7];

    double R13=(Rz[0]*Ry[0]+Rz[1]*Ry[3]+Rz[2]*Ry[6])*Rx[2]+(Rz[0]*Ry[1]+Rz[1]*Ry[4]+Rz[2]*Ry[7])*Rx[5]+(Rz[0]*Ry[2]+Rz[1]*Ry[5]+Rz[2]*Ry[8])*Rx[8];
    double R23=(Rz[3]*Ry[0]+Rz[4]*Ry[3]+Rz[5]*Ry[6])*Rx[2]+(Rz[3]*Ry[1]+Rz[4]*Ry[4]+Rz[5]*Ry[7])*Rx[5]+(Rz[3]*Ry[2]+Rz[4]*Ry[5]+Rz[5]*Ry[8])*Rx[8];
    double R33=(Rz[6]*Ry[0]+Rz[7]*Ry[3]+Rz[8]*Ry[6])*Rx[2]+(Rz[6]*Ry[1]+Rz[7]*Ry[4]+Rz[8]*Ry[7])*Rx[5]+(Rz[6]*Ry[2]+Rz[7]*Ry[5]+Rz[8]*Ry[8])*Rx[8];

        for(int i=start; i<(n_total)/3; ++i)
            {
                pos[3*i]  = pos[3*i]-o[0];
                pos[3*i+1]= pos[3*i+1]-o[1];
                pos[3*i+2]= pos[3*i+2]-o[2];

                double temp1 = R11*pos[3*i]+R12*pos[3*i+1]+R13*pos[3*i+2];
                double temp2 = R21*pos[3*i]+R22*pos[3*i+1]+R23*pos[3*i+2];
                double temp3 = R31*pos[3*i]+R32*pos[3*i+1]+R33*pos[3*i+2];
                pos[3*i] = temp1;
                pos[3*i+1] = temp2;
                pos[3*i+2] = temp3;

                pos[3*i]  = pos[3*i]+o[0];
                pos[3*i+1]= pos[3*i+1]+o[1];
                pos[3*i+2]= pos[3*i+2]+o[2];
            }


}


int interpBathymetry(double* sTrue, int* n, double xa, double xb, double ya, double yb, double height0, double hFreeSurface,
  int Nx, int Ny, double* bath, std::vector<double>& posFree, std::vector<double>& posFixed, double perturbation)
  {
    int nFreeTotal = 0;
    double dx = (xb-xa)/(double)Nx;
    double dy = (yb-ya)/(double)Ny;
    int i,j,k,l, k_low,l_low, k_up, l_up;
    double x, y, x_up, x_low, y_up, y_low;
    for(int i=0;i<=n[0];i++)
    {
      k = (int)floor( ( ((double)i) + 0.5)*sTrue[0]/dx);
      if(k< Nx)
      {
        k_low = k; k_up = k+1;
      }
      else if(k == Nx)
      {
        //should not appear as we stop s/2 before end of domain
        k_low = k-1; k_up = k;
      }
      else
      {
        printf("Error reading bathemetry\n");
      }
      x = xa + ((double)i+0.5)*sTrue[0];
      for(int j=0;j<=n[1];j++)
      {
        l = (int)floor( ( ((double)j)+0.5)*sTrue[1]/dy );
        if(l < Ny)
        {
          l_low = l; l_up = l+1;
        }
        else if(l == Ny)
        {
          //should not appear as we stop s/2 before end of domain
          l_low = l-1; l_up = l;
        }
        else{printf("Error reading bathemetry\n");
      }

      y = ya + ((double)j+0.5)*sTrue[1];

      //Bathemetric nodes around mesh nodes (x,y)
      x_up = xa + ((double)k_up)*dx;
      x_low = xa + ((double)k_low)*dx;
      y_up = ya + ((double)l_up)*dy;
      y_low = ya + ((double)l_low*dy);
      double z =( bath[k_low*(Ny+1) + l_low]*(x_up - x)*(y_up - y)
      + bath[k_up*(Ny+1) + l_low]*(x - x_low)*(y_up-y)
      + bath[k_low*(Ny+1) + l_up]*(x_up - x)*(y-y_low)
      + bath[k_up*(Ny+1) + l_up]*(x - x_low)*(y-y_low) )/(dx*dy) + height0;


      // generates number in the range -s*perturbation % and s*perturbation %
      std::default_random_engine generator; // A SEED MUST BE USE TO CHANGE VALUE AT EACH CALL
      std::uniform_real_distribution<double> distribution(-sTrue[0]*perturbation/100,sTrue[0]*perturbation/100);
      for(int m = 0; m <= n[2]; m++)
      {
        posFixed.push_back(x + distribution(generator));
        posFixed.push_back(y + distribution(generator));
        posFixed.push_back(z - ((double)m)*sTrue[2] + distribution(generator));
      }
      int nzFree = std::max((int) round( (hFreeSurface - z)/sTrue[2]),0); //round because it is allowed to go a bit up of the free surface !

      //Implement a possible difference for szFree and szFixed !!
      for(int p = 1; p <= nzFree;p++)
      {
        posFree.push_back(x + distribution(generator));
        posFree.push_back(y + distribution(generator));
        posFree.push_back(z + ((double)p)*sTrue[2] + distribution(generator));
      }
      nFreeTotal += nzFree;
    }
  }
  return nFreeTotal;
}

Error meshBathymetry(char* batFile, int numberGroundParticles, double height0, double hFreeSurface, double s, std::vector<double> &posFree,std::vector<double> &posFixed,  int* nPartFree, int* nPartFixed, double* volPart,
     double perturbation, bool stack)
     {

             double* bath;
             char buf[1000];
             // Opening files
             FILE *fp_bat;
             fp_bat = fopen(batFile, "r");

             if (fp_bat == NULL)
             {
               std::cout << "Bathymetry file not valid" << std::endl;
               return argumentError;
             }

             // Reading bathymetry parameters
             unsigned int bytesRead;
             bytesRead = fread(&buf, 8, 1, fp_bat);
             double xa = (*(double*)buf);
             //std::cout << xa << std::endl;
             bytesRead = fread(&buf, 8, 1, fp_bat);
             double xb = (*(double*)buf);
             //std::cout << xb << std::endl;
             bytesRead = fread(&buf, 8, 1, fp_bat);
             double ya = (*(double*)buf);
             //std::cout << ya << std::endl;
             bytesRead = fread(&buf, 8, 1, fp_bat);
             double yb = (*(double*)buf);
             //std::cout << yb << std::endl;
             bytesRead = fread(&buf, 4, 1, fp_bat);
             int Nx = (*(int*)buf);
             //std::cout << Nx << std::endl;
             bytesRead = fread(&buf, 4, 1, fp_bat);
             int Ny = (*(int*)buf);
             //std::cout << Ny << std::endl;


             // Reading the bathymetry
             bath = (double*) malloc((Nx+1)*(Ny+1)*sizeof(double));
             for(int i=0; i < ((Nx+1)*(Ny+1)); ++i)
             {
               bytesRead = fread(&buf, 8, 1, fp_bat);
               bath[i] = (*(double*)buf);
             }
             fclose(fp_bat);
             // Interpolating the bathymetry
             // Number of grid points -1 along x and y
             // Number of fixed particles along z
             int n[3];
             n[0] = (int) floor( (xb-xa)/s -(stack == true));
             n[1] = (int) floor( (yb-ya)/s -(stack == true));
             n[2] = numberGroundParticles -1 ;


             double sTrue[3];
             sTrue[0] = (xb-xa)/((double) (n[0]+(stack == true)));
             sTrue[1] = (yb-ya)/((double) (n[1]+(stack == true)));
             sTrue[2] = s;
             *volPart = sTrue[0]*sTrue[1]*sTrue[2];
             *nPartFixed = (n[0]+1)*(n[1]+1)*(n[2]+1);

             // memory allocation
             posFixed.reserve(posFixed.size() + *nPartFixed*3);
             // Impossible to know the posFree size at this point

             *nPartFree = interpBathymetry(sTrue, n, xa, xb, ya, yb,height0, hFreeSurface, Nx,Ny,bath, posFree, posFixed, perturbation);
             free(bath);

             return noError;
     }
