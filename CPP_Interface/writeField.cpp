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
    scalarsFree["pressureFree"] = &field->pressureFree;
    scalarsFree["densityFree"]  = &field->densityFree;
    vectorsFree["velocityFree"] = &field->speedFree;

    // Map Moving particules
    std::map<std::string, std::vector<double> *> scalarsMoving;
    std::map<std::string, std::vector<double> *> vectorsMoving;
    scalarsMoving["pressureMoving"] = &field->pressureMoving;
    scalarsMoving["densityMoving"]  = &field->densityMoving;
    vectorsMoving["velocityMoving"] = &field->speedMoving;

    // Map Fixed particules
    std::map<std::string, std::vector<double> *> scalarsFixed;
    std::map<std::string, std::vector<double> *> vectorsFixed;
    scalarsFixed["pressureFixed"] = &field->pressureFixed;
    scalarsFixed["densityFixed"]  = &field->densityFixed;

    // Save results to disk
    paraView("resultsFree", t, field->posFree, scalarsFree, vectorsFree);

}
