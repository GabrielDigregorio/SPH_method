#include "Main.h"
#include "Interface.h"

/*
 * In: field = stucture containing value to write
 *     t = time corresponding to the file to write
 * Out: speed_t.vtk, pos_t.vtk, ...
 */
void writeField(Field* field, double t)
{
    // Map Free particules
    std::map<std::string, std::vector<double> *> scalarsFree;
    std::map<std::string, std::vector<double> *> vectorsFree;
    scalars["pressureFree"] = &Field->pressureFree;
    scalars["densityFree"]  = &Field->densityFree;
    vectors["velocityFree"] = &Field->speedFree;

    // Map Moving particules
    std::map<std::string, std::vector<double> *> scalarsMoving;
    std::map<std::string, std::vector<double> *> vectorsMoving;
    scalars["pressureMoving"] = &Field->pressureMoving;
    scalars["densityMoving"]  = &Field->densityMoving;
    vectors["velocityMoving"] = &Field->velocityMoving;

    // Map Fixed particules
    std::map<std::string, std::vector<double> *> scalarsFixed;
    std::map<std::string, std::vector<double> *> vectorsFixed;
    scalars["pressureFixed"] = &Field->pressureFixed;
    scalars["densityFixed"]  = &Field->densityFixed;
    vectors["velocityFixed"] = &Field->velocityFixed;

    // Save results to disk
    paraview("resultsFree", t, posFree, scalarsFree, vectorsFree);

}
