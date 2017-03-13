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
  /// READING INPUT FILE
  std::string parameterFilename;
 std::string geometryFilename;
  std::string experimentFilename = "result"; // default name

  // Check the arguments
  if(argc<3) // Not enough argument
  {
    std::cout << "Invalid input files.\n";
    return EXIT_FAILURE;
  }
  else if (argc<4) // Use default name for the experiment (result)
    {parameterFilename = argv[1]; geometryFilename = argv[2];}
  else // Use default name for the experiment (result)
    {parameterFilename = argv[1]; geometryFilename = argv[2]; experimentFilename = argv[3];}

  //Read parameters
  Parameter* parameter =  new Parameter();
  readParameter(parameterFilename, parameter);

  //Read geometry
  Field* currentField =  new Field();
  readGeometry(geometryFilename, currentField);

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
  Field *nextField = new Field(), *tmpField = new Field();

  // Reserv Memory for each field of nextField and tmpField
  sizeField(nextField, currentField->nTotal);
  sizeField(tmpField, currentField->nTotal);

  unsigned int nMax = (unsigned int) ceil(parameter->T/parameter->k); //Validité de cette ligne à vérifier

  writeField(currentField, 0.0, parameter->format);
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
      writeField(nextField, n,  parameter->format);
      writeCount++;
    }
    tmpField = currentField;
    currentField = nextField;
    nextField = tmpField;
    std::cout << "----END time step #----" << n << "\n \n";
  }

  // Free all vectors and structurs
  cleanField(currentField);
  cleanField(nextField);
  cleanField(tmpField);
  cleanParameter(parameter);
  boxClear(boxes);
  boxClear(surrBoxesAll);

  return 0;
}

