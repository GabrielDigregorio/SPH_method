#include "Main.h"
#include "Physics.h"
#define M_PI           3.14159265358979323846  /* pi */
/*
 * In: field = structure containing the speed of particules (among others)
 *     parameter = structure containing the parameter usefull to know the movement of the wall
 * Out: Mise à jour des vitesses des parois mobiles
 */

 void updateMovingSpeed(Field* field, Parameter* parameter, double t)
 {
   // To be totally changed due to new strucure of the code !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Should compute speed particle by particle instead of for the whole domain !!!!
  /* int start;
   int end;
   start = field->nFixed + field->nFree;
   end = field->nTotal;

   switch(speedLaw)
   {
     case sine:
     {
       // uniform periodic movement in each direction
       int end=p+nx*ny*nz;
       for(int i=p ; i<end ; i++)
       {
         // x
         field->speed[0][i]=movx*sin(parameter->charactTime*t);
         // y
         field->speed[1][i]=movy*sin(parameter->charactTime*t);
         // z
         field->speed[2][i]=movz*sin(parameter->charactTime*t);
       }
     }
     break;
     case constant:
     {
       for(int i=start ; i<end ; i++)
       {
         // x
         field->speed[0][i]=movx;
         // y
         field->speed[1][i]=movy;
         // z
         field->speed[2][i]=movz;
       }
     }
     break;
     case exponential:
     std::cout << "Exponential moving not yet implemented\n"; // ATTENTION !!!
     break;
   }
   */
 }
