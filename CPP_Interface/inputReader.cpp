///**************************************************************************
/// SOURCE: Functions to read the parameters and geometry given by the user.
///**************************************************************************
#include "Physics.h"
#include "Main.h"
#include "Interface.h"

#include <sstream>
#include <algorithm>

#define N_UL 3
#define N_DATA 17
#define N_DATA_BATHYMETRY 5
#define N_PARAM 23

enum geomType{cube,cylinder,sphere};

/*
*Input:
*- type: Brick type to be read
*- inFile: Pointer to the input stream associated to the geometry file
*- currentField: Pointer to the structure to fill
*- posFree/posFree/posMoving: vector to store the position of the particles generated during brick reading
*- volVectorFree/Fixed/Moving: vector to store the volume of the particles generated during brick reading
*Decscription:
*Read a brick from the geometry file and generate the position and the volume of the particle inside the brick and store these information in the corresponding vectors.
*/

// Why "int type" and not "GoeomType type" ? would help identify the use of this argument ?
// Why passing vector with a pointer ? Isn't it possible with a reference ?
// Field* currentField useless here ?

Error readBrick(int type, std::ifstream* inFile, Parameter* parameter, std::vector<double>* posFree,
        std::vector<double>* posMoving, std::vector<double>* posFixed,
        std::vector<double>* volVectorFree,  std::vector<double>* volVectorFixed, std::vector<double>* volVectorMoving,std::vector<int>* typeFree,  std::vector<int>* typeFixed, std::vector<int>* typeMoving, int* numberMovingBoundaries)
{
    std::string buf;
    int cnt=0;
    char valueArray[1024];
    float brickData[N_DATA];

    while(cnt!=N_DATA && inFile->peek() != std::ifstream::traits_type::eof())
    {
        std::getline(*inFile, buf);
        if(1==sscanf(buf.c_str(),"%*[^#]#%s", valueArray)){
            std::cout << "Missing an element parameter.\n" << std::endl;
            return geometryError;
        }
        if(1==sscanf(buf.c_str(),"%*[^=]=%s", valueArray))
        {
            brickData[cnt]=atof(valueArray);
            ++cnt;
        }
        else{continue;}
    }
    if(cnt!=N_DATA)
    {
        std::cout << "Reached end of file before geometry reading.\n" << std::endl;
        return geometryError;
    }
    int c=(int)brickData[0];
    float s=brickData[1];
    float r=brickData[2];
    double o[3] = {brickData[3],brickData[4],brickData[5]};
    double L[3] = {brickData[6],brickData[7],brickData[8]};
    double teta[3]= {brickData[9],brickData[10],brickData[11]};
    int speedLaw= brickData[12];
    int charactTime = brickData[13];
    double movingDirection[3]={brickData[14],brickData[15],brickData[16]};
    for(int j = 0; j<3; j++)
    {
      if(L[j] < s)
      {
        std::cout << "Box smaller than s" << '\n';
        return geometryError;
      }
    }
    int nPart;
    double volPart;

    switch(c)
    {
        case freePart :
        switch(type)
        {
            case cube :
            meshcube(o, L,teta,s, *posFree, &nPart, &volPart, r, true);
            break;
            case cylinder :
            meshcylinder(o, L, s, *posFree, &nPart, &volPart, r, true);
            break;
            case sphere :
            meshsphere(o, L, s, *posFree, &nPart, &volPart, r, true);
            break;
        }
        for(cnt=0; cnt<nPart; ++cnt)
        {
            (*volVectorFree).push_back(volPart);
            (*typeFree).push_back(c);
        }
        break;
        case fixedPart :
        switch(type)
        {
            case cube :
            meshcube(o, L,teta, s, *posFixed, &nPart, &volPart, r, true);
            break;
            case cylinder :
            meshcylinder(o, L, s, *posFixed, &nPart, &volPart, r, true);
            break;
            case sphere :
            meshsphere(o, L, s, *posFixed, &nPart, &volPart, r, true);
            break;
        }
        for(cnt=0; cnt<nPart; ++cnt)
        {
            (*volVectorFixed).push_back(volPart);
            (*typeFixed).push_back(c);
        }
        break;
        case movingPart :
        {
        int IDMovingBoundary;
        *numberMovingBoundaries++;
        IDMovingBoundary = *numberMovingBoundaries - 1;
        
        for(int j=0 ; j<3 ; j++)
        {
            parameter->teta[j].push_back(teta[j]);
            parameter->movingDirection[j].push_back(movingDirection[j]);
        }
        parameter->charactTime.push_back(charactTime);
        parameter->speedLaw.push_back( (SpeedLaw) speedLaw);
        switch(type)
        {
            case cube :
            meshcube(o, L,teta,s, *posMoving, &nPart, &volPart, r, true);
            break;
            case cylinder :
            meshcylinder(o, L, s, *posMoving, &nPart, &volPart, r, true);
            break;
            case sphere :
            meshsphere(o, L, s, *posMoving, &nPart, &volPart, r, true);
            break;
        }
        for(cnt=0; cnt<nPart; ++cnt)
        {
            (*volVectorMoving).push_back(volPart);
            (*typeMoving).push_back(IDMovingBoundary + 2);  //Indeed, type = 0 is free, type = 1 is fixed and type > 1 is movingS !
        }
      }
        break;
    }
    return noError;
}

/*
*Input:
*- inFile: Pointer to the input stream associated to the geometry file
*- currentField: Pointer to the structure to fill
*- posFree/posFree/posMoving: vector to store the position of the particles generated during brick reading
*- volVectorFree/Fixed/Moving: vector to store the volume of the particles generated during brick reading
*Decscription:
*/

void readBathymetry(std::ifstream* inFile, std::vector<double>* posFree, std::vector<double>* posFixed,
        std::vector<double>* volVectorFree,  std::vector<double>* volVectorFixed, std::vector<int>* typeFree,  std::vector<int>* typeFixed)
{
    std::string buf;
    char batFile[64];
    int cnt=0;
    char valueArray[1024];
    float brickData[N_DATA_BATHYMETRY];

    std::getline(*inFile, buf);
    sscanf(buf.c_str(),"%*[^=]=%s", batFile);

    std::cout << "batFile name read: " << batFile << '\n';

    while(cnt!=N_DATA_BATHYMETRY)
    {
        std::getline(*inFile, buf);
        if(1==sscanf(buf.c_str(),"%*[^=]=%s", valueArray))
        {
            brickData[cnt]=atof(valueArray);
            ++cnt;
        }
        else{continue;}
    }

    float s=brickData[0];
    float r=brickData[1];
    double numberGroundParticles = (int) brickData[2];
    double height0 = brickData[3];
    double hFreeSurface = brickData[4];
    int nPartFree, nPartFixed;
    double volPart;

    std::cout << "BathyBrick read." << '\n';

    meshBathymetry(batFile, numberGroundParticles, height0, hFreeSurface, s, *posFree, *posFixed, &nPartFree, &nPartFixed, &volPart,
         r, true);

         std::cout << "Bathy meshed." << '\n';

         for(int i=0; i<nPartFree; i++)
         {
             (*volVectorFree).push_back(volPart);
             (*typeFree).push_back(freePart);
         }
         for(int i=0; i<nPartFixed; i++)
         {
             (*volVectorFixed).push_back(volPart);
             (*typeFixed).push_back(fixedPart);
         }
}

/*
*Input:
*- filename: name of the geometry file to read
*- volVector: vector to store the volume of the particles generated during reading the whole geometry
*Output:
*- 1: error
*- 0: no error
*Decscription:
*Read a entire geometry file and generate the position and the volume of the particles and store these informations in a structure.
*/
Error readGeometry(std::string filename, Field* currentField, Parameter* parameter, std::vector<double>* volVector)
{
    int numberMovingBoundaries = 0;
    std::vector<double> posFree, posFixed, posMoving, volVectorFree, volVectorFixed, volVectorMoving;
    std::vector<int> typeFree, typeFixed, typeMoving;
    std::ifstream inFile(filename);
    std::string buf;
    char valueArray[1024];
    int cnt=0;
    std::getline(inFile, buf);
    if(buf!="#GEOM1")
    {
        std::cout << "Wrong Geometry file type.\n" << std::endl;
        return geometryError;
    }
    while(inFile.peek() != std::ifstream::traits_type::eof())
    {
        std::getline(inFile, buf);
        buf.erase(std::remove(buf.begin(), buf.end(), ' '),buf.end());
        if(buf.size()!=0){
            switch(buf[0])
            {
                case '%' :
                continue;
                break;
                case '#' :
                buf.erase (0,1);
                if(buf=="domsz")
                {
                    while(cnt!=N_UL && inFile.peek() != std::ifstream::traits_type::eof())
                    {
                        std::getline(inFile, buf);
                        if(1==sscanf(buf.c_str(),"%*[^#]#%s", valueArray)){
                            std::cout << "Missing a domain parameter.\n" << std::endl;
                            return geometryError;
                        }
                        if(1==sscanf(buf.c_str(),"%*[^=]=%s", valueArray))
                        {
                            currentField->u[cnt]=atof(valueArray);
                            ++cnt;
                        }
                        else
                            continue;
                }
                if(cnt!=N_UL)
                {
                    std::cout << "Reached end of file before geometry reading.\n" << std::endl;
                    return geometryError;
                }
                cnt = 0;
                    while(cnt!=N_UL && inFile.peek() != std::ifstream::traits_type::eof())
                    {
                        std::getline(inFile, buf);
                        if(1==sscanf(buf.c_str(),"%*[^#]#%s", valueArray)){
                            std::cout << "Missing a domain parameter.\n" << std::endl;
                            return geometryError;
                        }
                        if(1==sscanf(buf.c_str(),"%*[^=]=%s", valueArray))
                        {
                            currentField->l[cnt]=atof(valueArray);
                            ++cnt;
                        }
                        else
                            continue;
                    }
                    if(cnt!=N_UL){
                        std::cout << "Reached end of file before geometry reading.\n" << std::endl;
                        return geometryError;
                    }
                }
                else if(buf=="brick")
                {
                    if(readBrick(cube,&inFile, parameter,
                        &posFree, &posMoving, &posFixed,
                        &volVectorFree, &volVectorFixed, &volVectorMoving,&typeFree, &typeFixed, &typeMoving, &numberMovingBoundaries)==geometryError){
                            return geometryError;
                    }
                }
                else if(buf=="cylin")
                {
                    if(readBrick(cylinder,&inFile, parameter,
                        &posFree, &posMoving, &posFixed,
                        &volVectorFree, &volVectorFixed, &volVectorMoving,&typeFree, &typeFixed, &typeMoving, &numberMovingBoundaries)==geometryError){
                            return geometryError;
                    }
                }
                else if(buf=="spher")
                {
                    if(readBrick(sphere,&inFile, parameter,
                        &posFree, &posMoving, &posFixed,
                        &volVectorFree, &volVectorFixed, &volVectorMoving,&typeFree, &typeFixed, &typeMoving, &numberMovingBoundaries)==geometryError){
                            return geometryError;
                    }
                }
                else if(buf=="bathy")
                {
                    std::cout << "Bathy recognized." << '\n';
                    readBathymetry(&inFile, &posFree, &posFixed,
                            &volVectorFree,  &volVectorFixed, &typeFree, &typeFixed);
                }
                else if(buf=="END_G")
                {
                    // Number of particles
                    currentField->nFree=posFree.size()/3;
                    currentField->nFixed=posFixed.size()/3;
                    currentField->nMoving=posMoving.size()/3;
                    currentField->nTotal= currentField->nFree + currentField->nFixed + currentField->nMoving;
                    // Position and volume vector sorting
                    posFree.insert(posFree.end(), posFixed.begin(), posFixed.end());
                    posFree.insert(posFree.end(), posMoving.begin(), posMoving.end());
                    volVectorFree.insert(volVectorFree.end(), volVectorFixed.begin(), volVectorFixed.end());
                    volVectorFree.insert(volVectorFree.end(), volVectorMoving.begin(), volVectorMoving.end());
                    typeFree.insert(typeFree.end(), typeFixed.begin(), typeFixed.end());
                    typeFree.insert(typeFree.end(), typeMoving.begin(), typeMoving.end());
                    (*volVector)=volVectorFree;

                    // Filling position vector
                    for(int i=0 ; i < currentField->nTotal ; i++)
                    {
                        currentField->type.push_back(typeFree[i]); //Possible to do this in a better way ?
                        for(int j=0 ; j<3 ; j++)
                        {
                            currentField->pos[j].push_back(posFree[3*i+j]);
                        }
                    }
                    return noError;
                }
                else
                {
                    std::cout <<"Unknown '"<<buf<<"' identifier.\n" << std::endl;
                    return geometryError;
                }
                break;
                default :
                    std::cout <<"Wrong tag in the geometry file.\n" << std::endl;
                    return geometryError;
                continue;
            }
        }
    }
    std::cout << "Reached end of file before geometry reading.\n" << std::endl;
    return geometryError;
}

/*
*Input:
*- filename: name of the parameter file to read
*- parameter: pointer the the structure to fill
*Output:
*- 1: error
*- 0: no error
*Decscription:
*Read a parameter file and the information in a structure.
*/
Error readParameter(std::string filename, Parameter* parameter)
{
    // open a file geometry.kzr
    std::ifstream inFile(filename);
    std::string buf;
    std::getline(inFile, buf);
    if(buf!="#PARA1")
    {
        std::cout << "Wrong parameter file type.\n" << std::endl;
        return parameterError;
    }
    while(inFile.peek() != std::ifstream::traits_type::eof())
    {
        std::getline(inFile, buf);
        buf.erase(std::remove(buf.begin(), buf.end(), ' '),buf.end());
        switch(buf[0])
        {
            case '%' :
            continue;
            break;
            case '#' :
            buf.erase (0,1);
            if(buf=="param")
            {
                int cnt=0;
                char valueArray[1024];
                while(cnt!=N_PARAM && inFile.peek() != std::ifstream::traits_type::eof())
                {
                    std::getline(inFile, buf);
                    if(1==sscanf(buf.c_str(),"%*[^=]=%s", valueArray))
                    {
                        if(cnt==0)
                            parameter->kh=atof(valueArray);
                        if(cnt==1)
                            parameter->k=atof(valueArray);
                        if(cnt==2)
                            parameter->T=atof(valueArray);
                        if(cnt==3)
                            parameter->densityRef=atof(valueArray);
                        if(cnt==4)
                            parameter->B=atof(valueArray);
                        if(cnt==5)
                            parameter->gamma=atof(valueArray);
                        if(cnt==6)
                            parameter->g=atof(valueArray);
                        if(cnt==7)
                            parameter->writeInterval=atof(valueArray);
                        if(cnt==8)
                            parameter->c=atof(valueArray);
                        if(cnt==9)
                            parameter->alpha=atof(valueArray);
                        if(cnt==10)
                            parameter->beta=atof(valueArray);
                        if(cnt==11)
                            parameter->epsilon=atof(valueArray);
                        if(cnt==12)
                            parameter->molarMass = atof(valueArray);
                        if (cnt==13)
                            parameter->temperature = atof(valueArray);
                        if (cnt==14)
                            parameter->theta = atof(valueArray);
                        if(cnt==15)
                        {
                          if( (0 <= atoi(valueArray)) && (atoi(valueArray) < NB_KERNEL_VALUE) )
                          {
                            parameter->kernel=(Kernel) atoi(valueArray);
                          }
                          else
                          {
                            std::cout <<"Invalid Kernel.\n" << std::endl;
                            return parameterError;
                          }
                        }
                        if(cnt==16)
                        {
                          if( (0 <= atoi(valueArray)) && (atoi(valueArray) < NB_VISCOSITY_VALUE) )
                          {
                            parameter->viscosityModel=(ViscosityModel) atoi(valueArray);
                          }
                          else
                          {
                            std::cout <<"Invalid viscosityModel.\n" << std::endl;
                            return parameterError;
                          }
                        }
                        if(cnt==17)
                        {
                          if( (0 <= atoi(valueArray)) && (atoi(valueArray) < NB_INTEGRATION_VALUE) )
                          {
                            parameter->integrationMethod=(IntegrationMethod) atoi(valueArray);
                          }
                          else
                          {
                            std::cout <<"Invalid integrationMethod.\n" << std::endl;
                            return parameterError;
                          }
                        }
                        if(cnt==18)
                        {
                          if( (0 <= atoi(valueArray)) && (atoi(valueArray) < NB_ADAPTATIVE_VALUE) )
                          {
                            parameter->adaptativeTimeStep=(AdaptativeTimeStep) atoi(valueArray);
                          }
                          else
                          {
                            std::cout <<"Invalid adaptativeTimeStep.\n" << std::endl;
                            return parameterError;
                          }
                        }
                        if(cnt==19)
                        {
                          if( (0 <= atoi(valueArray)) && (atoi(valueArray) < NB_DENSITYINIT_VALUE) )
                          {
                            parameter->densityInitMethod=(DensityInitMethod) atoi(valueArray);
                          }
                          else
                          {
                            std::cout <<"Invalid densityInitMethod.\n" << std::endl;
                            return parameterError;
                          }
                        }
                        if(cnt==20)
                        {
                          if( (0 <= atoi(valueArray)) && (atoi(valueArray) < NB_STATEEQUATION_VALUE) )
                          {
                            parameter->stateEquationMethod=(StateEquationMethod) atoi(valueArray);
                          }
                          else
                          {
                            std::cout <<"Invalid stateEquationMethod.\n" << std::endl;
                            return parameterError;
                          }
                        }
                        if(cnt==21)
                        {
                          if( (0 <= atoi(valueArray)) && (atoi(valueArray) < NB_MASSINIT_VALUE) )
                          {
                            parameter->massInitMethod=(MassInitMethod) atoi(valueArray);
                          }
                          else
                          {
                            std::cout <<"Invalid massInitMethod.\n" << std::endl;
                            return parameterError;
                          }
                        }
                        if(cnt==22)
                        {
                          if( (0 <= atoi(valueArray)) && (atoi(valueArray) < NB_FORMAT_VALUE) )
                          {
                            parameter->format=(Format) atoi(valueArray);
                          }
                          else
                          {
                            std::cout <<"Invalid outputFormat.\n" << std::endl;
                            return parameterError;
                          }
                        }

                        ++cnt;
                    }
                    else{continue;}
                }
                if(cnt!=N_PARAM)
                {
                    std::cout << "Reached end of file before parameters reading.\n" << std::endl;
                    return parameterError;
                }
            }
            else if(buf=="END_F")
            {
                return noError;
            }
            else
            {
                std::cout <<"Unknown '"<<buf<<"' identifier. \n" << std::endl;
                return parameterError;
            }
            break;
            default :
            continue;
        }
    }
}
