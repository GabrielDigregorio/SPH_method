#include "Main.h"
#include "Interface.h"
#include "Physics.h"
#include "Tools.h"

/*
* In: -argv[1]: path to the parameter file
-argv[2]: path to the geometry file
*/

int main(int argc, char *argv[])
{
  // READING INPUT FILE
  // Check argument file
  char* parameterFilename;
  char* geometryFilename;
  if(argc<3)
  {
    std::cout << "Invalid input files.\n";
    return EXIT_FAILURE;
  }
  else{parameterFilename = argv[1]; geometryFilename = argv[2];}

  //Read parameters
  Parameter* parameter;

  //To implement
  //readParameter(parameterFilename,parameter);

  //Read geometry
  Field* currentField;

  //To implement
  //readGeometry(geometryFilename,currentField);

  // INITIALISATION
  // SPEEDS
  currentField->speed.assign(currentField->nFree+currentField->nFixed+currentField->nMoving,0.0);
  //On pourrait aussi imaginer implementer une fonction speedInit(currentField,parameter) si on veut pouvoir partir d'un autre état que le repos, ici on initialise toutes les vitesses à zeros et seul les moving boundaries seront éventuellement modifiées


  // Initialisation des moving boundaries
  if(currentField->nMoving != 0){updateMovingSpeed(currentField,parameter,0.0);}

  // DENSITIES
  //To implement densityInit, use formula from Goffin p122
  densityInit(currentField, parameter);


  // PRESSURES
  //(might not be necessary to store them as they only depend on density?)(to discuss)
  //To implement use formala 3.39 from Goffin
  pressureComputation(currentField,parameter);

  //MASSES
  //To implement use formula 3.59 from Goffin
  massInit(currentField,parameter);

  // UPDATE & WRITTING
  Field *nextField, *tmpField;
  unsigned int nMax = (unsigned int) ceil(parameter->T/parameter->k); //Validité de cette ligne à vérifier
  //To implement, the value "0" stands for the time a which we write
  writeField(currentField, 0, Matlab);
  unsigned int writeCount = 1;

  bool reBoxing = true;

  // Creates the box mesh and describes their adjacent relations
  std::vector<std::vector<int> > boxes;
  std::vector<std::vector<int> > surrBoxesAll;

  for(unsigned int n = 1;n<=nMax;n++)
  {
    std::cout << "----BEGIN time step #----" << n << "\n \n";
    if(reBoxing == true)
    {
      boxMesh(currentField->l, currentField->u, parameter->kh, boxes, surrBoxesAll);
    }
    reBoxing = timeIntegration(currentField,nextField,parameter,boxes,surrBoxesAll,n);

    if(writeCount*parameter->writeInterval <= n*parameter->k)
    {
      writeField(currentField, n, Matlab);
      writeCount++;
    }
    tmpField = currentField;
    currentField = nextField;
    nextField = tmpField;
    std::cout << "----END time step #----" << n << "\n \n";
  }
  return 0;
}
