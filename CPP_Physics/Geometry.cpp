#include "Main.h"
#include "Physics.h"

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
        k_low = k-1; k_up = k;
        std::cout << "Problem!" << '\n';
      }
      else
      {
        printf("Error reading bathemetry\n");
      }
      x = xa + ((double)i+0.5)*sTrue[0];
      for(int j=0;j<n[1];j++)
      {
        l = (int)floor( ( ((double)j)+0.5)*sTrue[1]/dy );
        if(l < Ny)
        {
          l_low = l; l_up = l+1;
        }
        else if(l == Ny)
        {
          l_low = l-1; l_up = l;
          std::cout << "Problem!" << '\n';
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
      int nzFree = std::max((int) floor( (hFreeSurface - z)/sTrue[2]),0);

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

// Build a brick of regulary aligned particles with the center of mass at o(x,y,z).
//  - o[3]: corner of the cube with the lowest (x,y,z) values ??? Center ???
//  - L[3]: edge lengths along x,y and z
//  - s: particle spacing
//  - pertubation: percentage of perturbation in position of particles (equal 0 by default)

void meshcube(double o[3], double L[3], double s, std::vector<double> &pos, int* nPart, double* volPart,
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
Error meshBathymetry(char* batFile, int numberGroundParticles, double height0, double hFreeSurface, double s, std::vector<double> &posFree,std::vector<double> &posFixed,  int* nPartFree, int* nPartFixed, double* volPart,
     double perturbation, bool stack)
     {


       std::cout << "mesh Bathymetry entered." << '\n';

             double* bath;
             char buf[1000];
             // Opening files
             FILE *fp_bat;
             fp_bat = fopen(batFile, "r");
             /*
             if (access(batFile, 0)==0)
              {
               std::cout << "Bathymetry file not valid" << std::endl;
               return argumentError;
             }


             */
             // Reading bathymetry parameters
             unsigned int bytesRead;
             bytesRead = fread(&buf, 8, 1, fp_bat);
             double xa = (*(double*)buf);
             std::cout << xa << std::endl;
             bytesRead = fread(&buf, 8, 1, fp_bat);
             double xb = (*(double*)buf);
             std::cout << xb << std::endl;
             bytesRead = fread(&buf, 8, 1, fp_bat);
             double ya = (*(double*)buf);
             std::cout << ya << std::endl;
             bytesRead = fread(&buf, 8, 1, fp_bat);
             double yb = (*(double*)buf);
             std::cout << yb << std::endl;
             bytesRead = fread(&buf, 4, 1, fp_bat);
             int Nx = (*(int*)buf);
             std::cout << Nx << std::endl;
             bytesRead = fread(&buf, 4, 1, fp_bat);
             int Ny = (*(int*)buf);
             std::cout << Ny << std::endl;


             // Reading the bathymetry
             bath = (double*) malloc((Nx+1)*(Ny+1)*sizeof(double));
             for(int i=0; i < ((Nx+1)*(Ny+1)); ++i)
             {
               bytesRead = fread(&buf, 8, 1, fp_bat);
               bath[i] = (*(double*)buf);
               if(bath[i] > 10)
               {
                std::cout << bath[i] << '\n';
               }
             }
             fclose(fp_bat);


             std::cout << "Bathymetry read." << '\n';
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
