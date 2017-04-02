///**************************************************************************
/// SOURCE: Functions related to MPI.
///**************************************************************************
#include "Main.h"
#include "Interface.h"
#include "Physics.h"
#include "Tools.h"
#include "Structures.h"


void scatterField(Field* globalField, Field* currentField){

    // Checks if the number of processor is not too high
    // size < 3kh
    // if not: easy solution -> exit with error code (?)
    //         better solution -> handle the problem and do not use the unsless process (but inform the user)

    // Scatter properly the domain AND the particles, with their whole information
    // For that, a call to processUpdate could be done...


    currentField->nFree = globalField->nFree;
    currentField->nFixed = globalField->nFixed;
    currentField->nMoving = globalField->nMoving;
    currentField->nTotal = globalField->nTotal;
    for (int i = 0; i < 3; i++){
      currentField->l[i] = globalField->l[i];
      currentField->u[i] = globalField->u[i];
      currentField->pos[i] = globalField->pos[i];
      currentField->speed[i] = globalField->speed[i];
    }
    currentField->mass = globalField->mass;
    currentField->type = globalField->type;
    currentField->density = globalField->density;
    currentField->pressure = globalField->pressure;

}

void gatherField(Field* globalField, Field* currentField){
    // Gather all the current fields into the global Field (for output file writing for example)

    globalField->nFree = currentField->nFree;
    globalField->nFixed = currentField->nFixed;
    globalField->nMoving = currentField->nMoving;
    globalField->nTotal = currentField->nTotal;
    for (int i = 0; i < 3; i++){
      globalField->l[i] = currentField->l[i];
      globalField->u[i] = currentField->u[i];
      globalField->pos[i] = currentField->pos[i];
      globalField->speed[i] = currentField->speed[i];
    }
    globalField->mass = currentField->mass;
    globalField->type = currentField->type;
    globalField->density = currentField->density;
    globalField->pressure = currentField->pressure;


}

void processUpdate(Field* currentField){
    // Sorts the particles
    // Sends the particles that leave the domain
    // Receives the particles that enter the domain
    // Sorts the particles
    // Sends the edges to halos of surrounding domains
    // Receives the halos from edges of surrounding domains

}

void timeStepFinding(Field* currentField){
    // Find on each process the maximum acceptable time step

    // Reduce this information among all processes
}
