#include "Main.h"
#include "Physics.h"
#define M_PI           3.14159265358979323846  /* pi */
/*
 * In: field = structure containing the speed of particules (among others)
 *     parameter = structure containing the parameter usefull to know the movement of the wall
 * Out: Mise à jour des vitesses des parois mobiles
 */

 void updateMovingSpeed(Field* field, Parameter* parameter, double t,int IDmovingboundary,int i )//,int i)
 {
   
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
         field->speed[0][i]=movX*sin(charactTime*t);
         //
         field->speed[1][i]=movY*sin(charactTime*t);
         // z
         field->speed[2][i]=movZ*sin(charactTime*t);

     }
     break;
     case constant:
     {

         // x
         field->speed[0][i]=movX;
         // y
         field->speed[1][i]=movY;
         // z
         field->speed[2][i]=movZ;

     }
     break;
     case exponential:
     std::cout << "Exponential moving not yet implemented\n"; // ATTENTION !!!
     break;
   }

 }
