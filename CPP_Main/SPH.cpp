#include "Main.h"
#include "Interface.h"
#include "Physics.h"
#include "Tools.h"
#include "Structures.h"


/*
* In: -argv[1]: path to the parameter file
-argv[2]: path to the geometry file
*/
int main(int argc, char *argv[])
{
  std::cout << "----BEGIN argument checking----\n" << std::endl;

  std::string parameterFilename;
  std::string geometryFilename;
  std::string experimentFilename;

  // Check the arguments
  if(argc<3) // Not enough argument
  {
    std::cout << "\t Invalid input files.\n";
    return EXIT_FAILURE;
  }
  else if (argc<4) // Use default name for the experiment (result)
    {parameterFilename = argv[1]; geometryFilename = argv[2]; experimentFilename = "result";} // default name
  else // Use default name for the experiment (result)
    {parameterFilename = argv[1]; geometryFilename = argv[2]; experimentFilename = argv[3];}

  std::cout << "----END argument checking----\n" << std::endl;

  //Read parameters
  Parameter* parameter =  new Parameter();
  readParameter(parameterFilename, parameter);

  //Read geometry
  Field* currentField =  new Field();
  readGeometry(geometryFilename, currentField);


  sizeField(currentField, currentField->nTotal);

  // INITIALISATION (all the vectors should have the right size here!)
  // SPEEDS
  speedInit(currentField,parameter);

  // DENSITIES
  densityInit(currentField, parameter);

  // PRESSURES
  pressureComputation(currentField,parameter);

  //MASSES
  massInit(currentField,parameter);

  // UPDATE & WRITTING
  Field *nextField = new Field();
  Field *tmpField;

  // Reserv Memory for each field of nextField and copy const parameters of currentField
  sizeField(nextField, currentField->nTotal);
  for (int i = 0; i < 3; i++)
  {
    nextField->l[i] = currentField->l[i];
    nextField->u[i] = currentField->u[i];
  }
  nextField->nFree = currentField->nFree;
  nextField->nFixed = currentField->nFixed;
  nextField->nMoving = currentField->nMoving;
  nextField->nTotal = currentField->nTotal;

 for(int i=0 ; i<currentField->nTotal ; i++)
    nextField->mass[i] = currentField->mass[i];

  std::cout << "----BEGIN time step #0"<< std::endl;
  unsigned int nMax = (unsigned int) ceil(parameter->T/parameter->k); //Validité de cette ligne à vérifier

  std::cout << "\t Number of time steps = " << nMax << "\n" << std::endl;

  // Creat directory to store data
  //experimentFilename = creatDirectory(experimentFilename);
  writeField(currentField, 0.0, parameter, parameterFilename, geometryFilename, experimentFilename);
  unsigned int writeCount = 1;

  bool reBoxing = true;
  std::cout << "----END time step #0"<< std::endl;

  // Creates the box mesh and describes their adjacent relations
  std::vector<std::vector<int> > boxes;
  std::vector<std::vector<int> > surrBoxesAll;

  for(unsigned int n = 1;n<=nMax;n++)
  {
    std::cout << "----BEGIN time step #" << n << "----" <<"\n \n";
    if(reBoxing == true)
    {
      std::cout << "\t Reboxing...\n" << std::endl;
      boxes.resize(0);// VERY BAD
      surrBoxesAll.resize(0);// VERY BAD
      boxMesh(currentField->l, currentField->u, parameter->kh, boxes, surrBoxesAll);
    }
    //std::cout << "Mass: " << nextField->mass[0]<<"\n";
    reBoxing = timeIntegration(currentField,nextField,parameter,boxes,surrBoxesAll,n);

    if(writeCount*parameter->writeInterval <= n*parameter->k)
    {
      writeField(nextField, n, parameter, parameterFilename, geometryFilename, experimentFilename);
      writeCount++;
    }
    tmpField = currentField;
    currentField = nextField;
    nextField = tmpField;
    std::cout << "----END time step #" << n << "----" <<"\n \n";
  }

  // Free all vectors and structurs
  cleanField(currentField);
  cleanField(nextField);
  cleanParameter(parameter);

  return 0;
}
