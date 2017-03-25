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
	double duration;

	// Creates an error flag
	Error errorFlag = noError;

	// CPU time repartition
	std::clock_t start;
	std::vector<double> timeInfo(6, 0.0);

	start = std::clock();
	std::cout << "Argument checking..." << std::endl;
	std::string parameterFilename;
	std::string geometryFilename;
	std::string experimentFilename;

	// Check the arguments
	if (argc<3) // Not enough argument
	{
		std::cout << "Invalid input files.\n" << std::endl;
		errorFlag = argumentError;
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

	// Read parameters
	Parameter parameterInstance;
	Parameter* parameter = &parameterInstance;
	errorFlag = readParameter(parameterFilename, parameter);
	if (errorFlag != noError)
	{
		return errorFlag;
	}

	// Creation of temporary volume vector
	std::vector<double> volVector;

	// Read geometry
	Field currentFieldInstance;
	Field* currentField = &currentFieldInstance;
	errorFlag = readGeometry(geometryFilename, currentField, &volVector); //Why sending adress of volVector ? Wouldn't it be more understable to pass it by reference ?
	if (errorFlag != noError)
	{
		return errorFlag;
	}

	// Checking consistency of user datas
	errorFlag = consistency(parameter, currentField);
	if (errorFlag != noError)
	{
		return errorFlag;
	}

	std::cout << "Done.\n" << std::endl;

	// Initialisation
	speedInit(currentField, parameter);
	densityInit(currentField, parameter);
	pressureInit(currentField, parameter);
	massInit(currentField, parameter, volVector);
	volVector.clear();
	
	// Creates field to store result of update
	Field nextFieldInstance;
	Field* nextField = &nextFieldInstance;
	copyField(currentField, nextField);

	unsigned int nMax = (unsigned int)ceil(parameter->T / parameter->k);
	std::cout << "Number of time steps = " << nMax << "\n" << std::endl;
	std::cout << "Number of free particles = " << currentField->nFree << "\n" << std::endl;
	std::cout << "Number of fixed particles = " << currentField->nFixed << "\n" << std::endl;
	std::cout << "Number of particles with imposed speed = " << currentField->nMoving << "\n" << std::endl;


	// Writes the initial configuration
	writeField(currentField, 0.0, parameter, parameterFilename, geometryFilename, experimentFilename);
	unsigned int writeCount = 1;

	// To mesh at least at the first time step
	bool reBoxing = true;

	// Creates the box mesh and describes their adjacent relations
	std::vector<std::vector<int> > boxes;
	std::vector<std::vector<int> > surrBoxesAll;

	timeInfo[0] = (std::clock() - start) / (double)CLOCKS_PER_SEC;

	// Loop on time
	std::cout << "Time integration progress:\n" << std::endl;
	std::cout << "0%----------------------------------------------100%\n[";

	for (unsigned int n = 1; n <= nMax; n++)
	{
		// Rebox the domain if h has sufficiently changed
		if (reBoxing == true)
		{
			start = std::clock();
			boxes.resize(0);// VERY BAD, TO CHANGE !!! How to do this properly ?
			surrBoxesAll.resize(0);// VERY BAD, TO CHANGE !!! How to do this properly ?
			boxMesh(currentField->l, currentField->u, parameter->kh, boxes, surrBoxesAll);
			timeInfo[1] += (std::clock() - start) / (double)CLOCKS_PER_SEC;
		}

		// Time integration
		reBoxing = timeIntegration(currentField, nextField, parameter, boxes, surrBoxesAll, n, timeInfo);

		// Write field when needed
		start = std::clock();
		if (writeCount*parameter->writeInterval <= n*parameter->k)
		{
			writeField(nextField, n, parameter, parameterFilename, geometryFilename, experimentFilename);
			writeCount++;
		}
		timeInfo[5] += (std::clock() - start) / (double)CLOCKS_PER_SEC;

		start = std::clock();

		// Swap two fields
		swapField(&currentField, &nextField);

		timeInfo[4] += (std::clock() - start) / (double)CLOCKS_PER_SEC;

		// Fancy progress bar
		if (!((50 * n) % nMax))
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
	std::cout << "\t- TOTAL  \t" << (std::clock() - startExperimentTimeClock) / (double)CLOCKS_PER_SEC << "\n";
	std::cout << "NB : Total - sum of times = time capture duration (!!)\n";

	return 0;
}
