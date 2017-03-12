#include "Main.h"
#include "Interface.h"
#include "Physics.h"
#include "Tools.h"
#include "Structures.h"

int main(int argc, char *argv[]){
        char* parameterFilename;
        char* geometryFilename;
        if(argc<3)
        {
                std::cout << "Invalid input files.\n";
            return EXIT_FAILURE;
        }
        else{parameterFilename = argv[1]; geometryFilename = argv[2];}


        //Read geometry
        Field currentField;
        readGeometry(geometryFilename, &currentField);

        //Parameter parameter;
        Parameter parameter;
        readParameter(parameterFilename, &parameter);
        return 0;
}
