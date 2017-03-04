#include "../Headers/SPH.hpp"
#include "../Headers/Playground.hpp"

// Structure qui doit être remplie lors de la lecture du fichier de paramètre (il faudra surement changer la place de cette déclaration)
/*
 * s = space interval
 * kh = smothing length
 * k = time step
 * T = simulation time
 * densityRef = density of the fluid at atmospheric pressure
 * l & u = lower and upper limit of the domain
 * B & gamma = fluid constants
 * g = gravity
 * writeInteval = time interval between two outputs file are generated
 * integrationMethod = euler ou RK2
 * densityInitMethod = hydrosatic, etc.
 * stateEquationMethod = quasiIncompressible, perfectGas, etc.
 * massInitMethod = violeau2012 (all particles have same volumes), etc.
 * speedLaw = To be determined, will dictate the behaviour of moving boundaries
*/
struct Parameter {
    double s, kh, k, T, densityRef, B, gamma, g, writeInterval;
    double l[3];
    double u[3];
    std::string integrationMethod, densityInitMethod, stateEquationMethod, massInitMethod, speedLaw;
};

// Structure qui doit être remplie lors de la lecture du fichier de géométrie (il faudra surement changer la place de cette déclaration aussi). Cette structure contient toute l'information utile de nos simulations.
struct Field {
    std::vector<double> posFree;
    std::vector<double> posMoving;
    std::vector<double> posFixed;

    std::vector<double> speedFree;
    std::vector<double> speedMoving;
    //Speed fixed = 0 of course

    std::vector<double> densityFree;
    std::vector<double> densityMoving;
    std::vector<double> densityFixed;

    std::vector<double> pressureFree;
    std::vector<double> pressureMoving;
    std::vector<double> pressureFixed;

    std::vector<double> massFree;
    std::vector<double> massMoving;
    std::vector<double> massFixed;
};

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
        readParameter(parameterFilename,parameter);

        //Read geometry
        Field* currentField;
        //To implement
        readGeometry(geometryFilename,currentField);

    // INITIALISATION
        // SPEEDS
            currentField->speedFree.assign(currentField->posFree.size(),0.0);
            //On pourrait aussi imaginer implementer une fonction speedInit(currentField,parameter) si on veut pouvoir partir d'un autre état que le repos


            //To implement (not to do in a first time because it is much complicated and we must discuss this.
            // the value 0 correpsonds to time t=0, the function will be reused for later times
            // On pourrait ne passer que parameter->speedLaw mais le passage par pointeur est tout aussi efficace et on a accès à tous les parametres dans le cas ou on en aurait besoin
            updateMovingSpeed(currentField,parameter,0.0);

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
            Field* nextField;
            unsigned int nMax = (unsigned int) ceil(T/k); //Validité de cette ligne à vérifier
            //To implement, the value "0" stands for the time a which we write
            writeField(currentField,0);
            unsigned int writeCount = 1;

            for(int n = 1;n<=nMax;n++)
            {
                timeIntegration(currentField,nextField,parameter,n);

                if(writeCount*parameter->writeInterval <= n*parameter->k)
                {
                    writeField(currentField,n*parameter->k);
                    writeCount++;
                }
            }
    return 0;
}
