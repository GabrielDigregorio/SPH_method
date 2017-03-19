#include "Main.h"
#include "Interface.h"
#include "Physics.h"
#include "Tools.h"
#include "Structures.h"

std::clock_t startExperimentTimeClock;

/*
* In: -argv[1]: path to the parameter file
-argv[2]: path to the geometry file
*/
int main(int argc, char *argv[]){
    // Record algorithm performance
    startExperimentTimeClock = std::clock();
    double duration;

    // CPU time repartition
    std::clock_t start;
    std::vector<double> timeInfo(6,0.0);

    start = std::clock();
    std::cout << "----BEGIN argument checking----\n" << std::endl;
    std::string parameterFilename;
    std::string geometryFilename;
    std::string experimentFilename;

    // Check the arguments
    if(argc<3){ // Not enough argument
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
    // Initialisation volume vector
    std::vector<double> volVector;
    //Read geometry
    Field* currentField =  new Field();
    readGeometry(geometryFilename, currentField, &volVector);
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
    for (int i = 0; i < 3; i++){
        nextField->l[i] = currentField->l[i];
        nextField->u[i] = currentField->u[i];
    }
    nextField->nFree = currentField->nFree;
    nextField->nFixed = currentField->nFixed;
    nextField->nMoving = currentField->nMoving;
    nextField->nTotal = currentField->nTotal;
    for(int i=0 ; i<currentField->nTotal ; i++)
        nextField->mass[i] = currentField->mass[i];
    unsigned int nMax = (unsigned int) ceil(parameter->T/parameter->k);
    std::cout << "\t Number of time steps = " << nMax << "\n" << std::endl;

    // Writes the initial configuration
    writeField(currentField, 0.0, parameter, parameterFilename, geometryFilename, experimentFilename);
    unsigned int writeCount = 1;

    // To mesh at least at the first time step
    bool reBoxing = true;

    // Creates the box mesh and describes their adjacent relations
    std::vector<std::vector<int> > boxes;
    std::vector<std::vector<int> > surrBoxesAll;

    timeInfo[0] = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

    // Loop on time
    std::cout << "----BEGIN time integration----\n"<< std::endl;
    std::cout << "0%----------------------------------------------100%\n[";

    for(unsigned int n = 1;n<=nMax;n++){
        // Rebox the domain if h has sufficiently changed
        if(reBoxing == true){
            start = std::clock();
            boxes.resize(0);// VERY BAD, TO CHANGE !!!
            surrBoxesAll.resize(0);// VERY BAD, TO CHANGE !!!
            boxMesh(currentField->l, currentField->u, parameter->kh, boxes, surrBoxesAll);
            timeInfo[1] += ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
        }
        // Time integration
        reBoxing = timeIntegration(currentField, nextField, parameter, boxes, surrBoxesAll, n, timeInfo);
        // Write field when needed
        start = std::clock();
        if(writeCount*parameter->writeInterval <= n*parameter->k){
            writeField(nextField, n, parameter, parameterFilename, geometryFilename, experimentFilename);
            writeCount++;
        }
        timeInfo[5] += ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

        start = std::clock();
        tmpField = currentField;
        currentField = nextField;
        nextField = tmpField;
        timeInfo[4] += ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

        // Fancy progress bar
        if(!((50*n)%nMax))
            std::cout << ">" << std::flush;
    }
    std::cout << "]\n" << std::endl;

    std::cout << "TIME INFORMATION:\n";
    std::cout << "\t- Initial\t" << timeInfo[0] << "\n";
    std::cout << "\t- Neighbors\t" << timeInfo[1] << "\n";
    std::cout << "\t- Continuity\t" << timeInfo[2] << "\n";
    std::cout << "\t- Momentum\t" << timeInfo[3] << "\n";
    std::cout << "\t- Update\t" << timeInfo[4] << "\n";
    std::cout << "\t- Writing\t" << timeInfo[5] << "\n";
    std::cout << "\t- TOTAL  \t" << ( std::clock() - startExperimentTimeClock ) / (double) CLOCKS_PER_SEC<<"\n";
    std::cout << "NB : Total - sum of times = time capture duration (!!)\n";

    // Free all vectors and structurs
    cleanField(currentField);
    cleanField(nextField);
    cleanParameter(parameter);

    return 0;
}
