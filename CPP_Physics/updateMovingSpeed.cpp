#include "Main.h"
#include "Physics.h"
#define M_PI           3.14159265358979323846  /* pi */
/*
 * In: field = structure containing the speed of particules (among others)
 *     parameter = structure containing the parameter usefull to know the movement of the wall
 * Out: Mise à jour des vitesses des parois mobiles
 */

 void updateMovingSpeed(Field* field, Parameter* parameter, double t,int indexI,int ii )
 {
   int IDmovingboundary=0;
    for(int i=0; i<parameter->Index.size();i++)
    {
      if(indexI>=parameter->Index[i])
      {
        IDmovingboundary=i;
      }
    }
    double movX=parameter->movingDirection[0][IDmovingboundary];
    double movY=parameter->movingDirection[1][IDmovingboundary];
    double movZ=parameter->movingDirection[2][IDmovingboundary];
    //std::vector<double> teta[3];
    double charactTime=parameter->charactTime[IDmovingboundary];
   switch(parameter->speedLaw[IDmovingboundary])
   {
     case sine:
     {
       // uniform periodic movement in each direction
         // x
         field->speed[0][ii]=movX*sin(charactTime*t);
         // 
         field->speed[1][ii]=movY*sin(charactTime*t);
         // z
         field->speed[2][ii]=movZ*sin(charactTime*t);
     }
     break;
     case constant:
     {    
         // x
         field->speed[0][ii]=movX;
         // y
         field->speed[1][ii]=movY;
         // z
         field->speed[2][ii]=movZ; 
     }
     break;
     case exponential:
     std::cout << "Exponential moving not yet implemented\n"; // ATTENTION !!!
     break;
     case level_arm:
     {  
       double tetax=(parameter->teta[0][IDmovingboundary]);
       double tetay=(parameter->teta[1][IDmovingboundary]);
       double tetaz=(parameter->teta[2][IDmovingboundary]);
       
       int identifier=indexI-2;
       int counter=0;
       while(counter<IDmovingboundary)
       {
        double s=parameter->spacingS[counter];
        double L1= parameter->Dimension[0][counter];
        double L2= parameter->Dimension[1][counter];
        double L3= parameter->Dimension[2][counter];
        int ni = int(ceil(L1/s));
        int nj = int(ceil(L2/s));
        int nk = int(ceil(L3/s));
        identifier=identifier-(ni*nk*nj);
        counter++;
       }
        double s=parameter->spacingS[IDmovingboundary];
        double L1= parameter->Dimension[0][IDmovingboundary];
        double L2= parameter->Dimension[1][IDmovingboundary];
        double L3= parameter->Dimension[2][IDmovingboundary];
        
        // calculate nb of particles along each direction from target size "s"
        int ni = int(ceil(L1/s));
        int nj = int(ceil(L2/s));
        int nk = int(ceil(L3/s));

        int levierI=identifier/(ni*nj);// give the lever arm , based on the construction block
        double levier=levierI/(double)20;
        double Ampli=parameter->ampliRota[IDmovingboundary];

        double angle=Ampli/180.0*M_PI*cos(charactTime*t);
        double angleDot=-charactTime*Ampli/180.0*M_PI*sin(charactTime*t);
        // in the basis relative to the block
        double vxRelBasis=levier*cos(angle)*angleDot;
        double vyRelBasis=0.0;
        double vzRelBasis=-levier*sin(angle)*angleDot;
        std::vector<double> Rx(9,0);
        std::vector<double> Ry(9,0);
        std::vector<double> Rz(9,0);
       // in hte real basis 
        Rx[0]=1;
        Rx[1]=0;
        Rx[2]=0;
        Rx[3]=0;
        Rx[4]=cos(tetax/180.0*M_PI);
        Rx[5]=-sin(tetax/180.0*M_PI);
        Rx[6]=0;
        Rx[7]=sin(tetax/180.0*M_PI);
        Rx[8]=cos(tetax/180.0*M_PI);
        
       /* Ry[0]=cos(tetay/180.0*M_PI);
        Ry[1]=0;
        Ry[2]=sin(tetay/180.0*M_PI);
        Ry[3]=0;
        Ry[4]=1;
        Ry[5]=0;
        Ry[6]=-sin(tetay/180.0*M_PI);
        Ry[7]=0;
        Ry[8]=cos(tetay/180.0*M_PI);*/
        Ry[0]=1;
        Ry[1]=0;
        Ry[2]=0;
        Ry[3]=0;
        Ry[4]=1;
        Ry[5]=0;
        Ry[6]=0;
        Ry[7]=0;
        Ry[8]=1;

        Rz[0]=cos(tetaz/180.0*M_PI);
        Rz[1]=-sin(tetaz/180.0*M_PI);
        Rz[2]=0;
        Rz[3]=sin(tetaz/180.0*M_PI);
        Rz[4]=cos(tetaz/180.0*M_PI);
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
        
        double temp1=R11*vxRelBasis+R13*vzRelBasis;
        double temp2=R21*vxRelBasis+R23*vzRelBasis;
        double temp3=R31*vxRelBasis+R33*vzRelBasis;
         // x
         field->speed[0][ii]=temp1;
         // y
         field->speed[1][ii]=temp2;
         // z
         field->speed[2][ii]=temp3;
        
     }
     break;
   }
   
 }

