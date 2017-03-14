#include "Main.h"
#include "Physics.h"

/*
 * In: field = structure containing the speed of particules (among others)
 *     parameter = structure containing the parameter usefull to know the movement of the wall
 * Out: Mise à jour des vitesses des parois mobiles
 */

void updateMovingSpeed(Field* field, Parameter* parameter, double t)
{
    int start;
    int end;

    switch(parameter->speedLaw){
    case sine:
        // uniform periodic movement in each direction
        start = field->nFixed + field->nFree;
        end = field->nTotal;

        for(int i=start ; i<end ; i++){
            // x
            field->speed[3*i]=parameter->movingDirection[0]*sin(parameter->charactTime*t);
            // y
            field->speed[3*i+1]=parameter->movingDirection[1]*sin(parameter->charactTime*t);
            // z
            field->speed[3*i+2]=parameter->movingDirection[2]*sin(parameter->charactTime*t);
        }
        break;
    case constant:
        std::cout << "Constant moving not yet implemented\n"; // ATTENTION !!!
        break;
    case exponential:
        std::cout << "Exponential moving not yet implemented\n"; // ATTENTION !!!
        break;
    }
}
