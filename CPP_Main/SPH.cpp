#include "Main.h"
#include "Interface.h"
#include "Physics.h"
#include "Tools.h"
#include <ctime>

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
            unsigned int nMax = (unsigned int) ceil(parameter->T/parameter->k); //Validité de cette ligne à vérifier
            //To implement, the value "0" stands for the time a which we write
            writeField(currentField,0);
            unsigned int writeCount = 1;

            for(unsigned int n = 1;n<=nMax;n++)
            {
                //timeIntegration(currentField,nextField,parameter,n);

                if(writeCount*parameter->writeInterval <= n*parameter->k)
                {
                    writeField(currentField,n*parameter->k);
                    writeCount++;
                }
            }
    return 0;
}
