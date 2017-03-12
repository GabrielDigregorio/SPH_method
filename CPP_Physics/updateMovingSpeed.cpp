#include "Main.h"
#include "Physics.h"

/*
 * In: field = structure containing the speed of particules (among others)
 *     parameter = structure containing the parameter usefull to know the movement of the wall
 * Out: Mise à jour des vitesses des parois mobiles
 */
void updateMovingSpeed(Field* field, Parameter* parameter, double t, MoveMod myMod)
{
    switch(myMod){

    case 1 :
        // uniforme periodic mouvement in each direction
        // case 1: simple translation right to left ( x direction ) ( sin pour tout )
        // depending on the parameter choose
        int nb_free;
        int nb_fixed;
        int nb_mov;
        nb_free=field->nFree;
        nb_fixed=field->nFixed;
        nb_mov=field->nMoving;
        int start=nb_fixed+nb_free;
        int end =nb_fixed+nb_free+nb_mov;
        // maybe put this in a vector of parameter ?
        double A1=1;
        double A2=1;
        double A3=1;
        double w1=1;
        double w2=1;
        double w3=1;
        for(int i=start ; i<end ; i=i+3)
        {
            // x
            field->speed[i]=A1*sin(w1*t);
            // y
            field->speed[i+1]=A2*sin(w2*t);
            // z
            field->speed[i+2]=A3*sin(w3*t);
        }

    break;

    case 2 : // case 2: mouvement levier avec le point le plus bas fixe.

    break;

    default :
    }
}
