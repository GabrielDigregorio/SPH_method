///**************************************************************************
/// SOURCE: Main function of the smoothed particle hydrodynamics (SPH) solver.
///**************************************************************************
#include "Main.h"
#include "Interface.h"
#include "Physics.h"
#include "Tools.h"
#include "Structures.h"

std::clock_t startExperimentTimeClock;

/*
*Input:
*- argv[1]: name of the input parameter file (mandatory)
*- argv[2]: name of the input geometry file (mandatory)
*- argv[3]: name of the output result (optional, default name is "result.txt")
*
*Decscription:
*Run the SPH solver for a given geometry and a given set of parameter and write the result in an output file.
*/
int main(int argc, char *argv[])
{
	// Record algorithm performance
	startExperimentTimeClock = std::clock();

	// MPI Initialization
	MPI_Init(&argc, &argv);
	int nTasks, procID;
	MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &procID);

	// Creates an error flag
	Error errorFlag = noError;

	// CPU time repartition information variables
	std::clock_t start;
	std::vector<double> timeInfo(6, 0.0);

	// Checks and gets the arguments
	if(procID==0){std::cout << "Initialization..." << std::endl;}
	std::string parameterFilename;
	std::string geometryFilename;
	std::string experimentFilename;
	if (argc<3) // Not enough arguments
	{
		std::cout << "Invalid input files.\n" << std::endl;
		errorFlag = argumentError;
		MPI_Finalize();
		return errorFlag;
	}
	else if (argc<4) // Use default name for the experiment (result)
	{
		parameterFilename = argv[1];
		geometryFilename = argv[2];
		experimentFilename = "result";
	}
	else {
		parameterFilename = argv[1];
		geometryFilename = argv[2];
		experimentFilename = argv[3];
	}

	// Main variables declaration
	Parameter parameterInstance;
	Parameter* parameter = &parameterInstance;
	Field currentFieldInstance;
	Field* currentField = &currentFieldInstance;
	Field nextFieldInstance;
	Field* nextField = &nextFieldInstance;
	Field globalFieldInstance; // Used by node 0 only
	Field* globalField = &globalFieldInstance; // Used by node 0 only

	start = std::clock();
	// Reads parameters (each process) and geometry (process 0) and checks their consistency
	errorFlag = readParameter(parameterFilename, parameter);
	if(errorFlag != noError){MPI_Finalize(); return errorFlag;}
	if(procID==0){
		errorFlag = initializeField(geometryFilename, globalField, parameter);
	}
	MPI_Bcast(&errorFlag, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if(errorFlag != noError){MPI_Finalize(); return errorFlag;}

	// Writes the initial configuration
	if(procID==0){writeField(globalField, 0.0, parameter, parameterFilename, geometryFilename, experimentFilename);}
	unsigned int writeCount = 1;

	// Scatters the globalField from node 0 into the currentField of all nodes
	scatterField(globalField, currentField);

	// Copies the invariant information about the field
	copyField(currentField, nextField);

	// Initialization done
	if(procID==0){std::cout << "Done.\n"  << std::endl;}

	// Information on the simulation
	if(procID==0){
		if (parameter->adaptativeTimeStep == no){
			unsigned int nMax = (unsigned int)ceil(parameter->T / parameter->k);
			std::cout << "Number of time steps = " << nMax << "\n" << std::endl;
		}
		else
	        std::cout << "Number of time steps = " << "not defined (adaptative time step)" << "\n" << std::endl;
		std::cout << "Number of free particles = " << currentField->nFree << "\n" << std::endl;
		std::cout << "Number of fixed particles = " << currentField->nFixed << "\n" << std::endl;
		std::cout << "Number of particles with imposed speed = " << currentField->nMoving << "\n" << std::endl;
	}

	// Declares the box mesh and their adjacent relations variables
	bool reBoxing = true; // To mesh at least at the first time step
	std::vector<std::vector<int> > boxes;
	std::vector<std::vector<int> > surrBoxesAll;

	timeInfo[0] = (std::clock() - start) / (double)CLOCKS_PER_SEC;

	// ------------ LOOP ON TIME ------------
	if(procID==0){
		std::cout << "Time integration progress:\n" << std::endl;
		std::cout << "0%-----------------------------------------------100%\n[";
	}
	unsigned int loadingBar = 0;
	double currentTime = 0.0; // Current time of the simulation
	for (unsigned int n = 1; currentTime < parameter->T; n++)
	{
		// Major MPI communication: the time integration is preparated
		processUpdate(currentField);

		currentField->nextK = parameter->k;

		// Rebox the domain if h has sufficiently changed (PROBABLY NOT TO BE LEFT HERE)
		if (reBoxing == true){
			start = std::clock();
			boxes.resize(0);// VERY BAD, TO CHANGE !!! How to do this properly ?
			surrBoxesAll.resize(0);// VERY BAD, TO CHANGE !!! How to do this properly ?
			boxMesh(currentField->l, currentField->u, parameter->kh, boxes, surrBoxesAll);
			timeInfo[1] += (std::clock() - start) / (double)CLOCKS_PER_SEC;
        }
		// Solve the time step
        reBoxing = timeIntegration(currentField, nextField, parameter, boxes, surrBoxesAll, currentTime,parameter->k, timeInfo);

		// Adaptative time step
		timeStepFinding(currentField);
		currentTime += parameter->k; // Temporary
		parameter->k = currentField->nextK; // Temporary

		// Swap two fields
		swapField(&currentField, &nextField);

		// Write field when needed
		start = std::clock();
    	if (writeCount*parameter->writeInterval <= currentTime)
		{
			gatherField(globalField, currentField);
			if(procID==0){writeField(currentField, n, parameter, parameterFilename, geometryFilename, experimentFilename);}
			writeCount++;
		}
		timeInfo[5] += (std::clock() - start) / (double)CLOCKS_PER_SEC;

		// Fancy progress bar (ok if at least 50 time step)
		if ( procID==0 && currentTime > loadingBar * parameter->T/50.0){
			std::cout << ">" << std::flush;
			loadingBar++;
		}
	}

	// MPI Finalize
	MPI_Finalize();

	if(procID==0){
		std::cout << "]\n" << std::endl;
		std::cout << "TIME INFORMATION:\n";
		std::cout << "\t- Initial\t" << timeInfo[0] << "\n";
		std::cout << "\t- Neighbors\t" << timeInfo[1] << "\n";
		std::cout << "\t- Continuity\t" << timeInfo[2] << "\n";
		std::cout << "\t- Momentum\t" << timeInfo[3] << "\n";
		std::cout << "\t- Update\t" << timeInfo[4] << "\n";
		std::cout << "\t- Writing\t" << timeInfo[5] << "\n";
		std::cout << "\t- TOTAL  \t" << (std::clock() - startExperimentTimeClock) / (double)CLOCKS_PER_SEC << "\n";
		std::cout << "NB : Total - sum of times = time capture duration (!!)\n";
	}

	return 0;
}
