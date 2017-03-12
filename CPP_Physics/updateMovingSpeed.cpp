#include "Main.h"
#include "Physics.h"

/*
 * In: field = structure containing the speed of particules (among others)
 *     parameter = structure containing the parameter usefull to know the movement of the wall
 * Out: Mise à jour des vitesses des parois mobiles
 */

void updateMovingSpeed(Field* field, Parameter* parameter, double t)
{
    switch(parameter->speedLaw){

    case sine :
        // uniforme periodic mouvement in each direction

        int start = field->nFixed + field->nFree;
        int end = field->nTotal;

        for(int i=start ; i<end-2 ; i=i+3)
        {
            // x
            field->speed[i]=parameter->movingDirection[1]*sin(parameter->charactTime*t);
            // y
            field->speed[i+1]=parameter->movingDirection[2]*sin(parameter->charactTime*t);
            // z
            field->speed[i+2]=parameter->movingDirection[3]*sin(parameter->charactTime*t);
        }

    break;

    //case 2 :

    //break;

    //default :
    }
}
