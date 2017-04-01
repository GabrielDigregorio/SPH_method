///**************************************************************************
/// SOURCE: Functions to read the parameters and geometry given by the user.
///**************************************************************************
#include "Physics.h"
#include "Main.h"
#include "Interface.h"

#include <sstream>
#include <algorithm>

#define N_UL 3
#define N_DATA 16
#define N_PARAM 24

enum geomType{cube,cylinder,sphere};
enum boundCondition{freePart, movingPart, fixedPart};

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

Error readBrick(int type, std::ifstream* inFile, Field* currentField, std::vector<double>* posFree,
        std::vector<double>* posMoving, std::vector<double>* posFixed,
        std::vector<double>* volVectorFree,  std::vector<double>* volVectorFixed, std::vector<double>* volVectorMoving)
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
    if(cnt!=N_DATA){
        std::cout << "Reached end of file before geometry reading.\n" << std::endl;
        return geometryError;
    }
    int c=(int)brickData[0];
    float s=brickData[1];
    float r=brickData[2];
    double o[3] = {brickData[3],brickData[4],brickData[5]};
    double L[3] = {brickData[6],brickData[7],brickData[8]};
    double tempo[3]= {brickData[6],brickData[7],brickData[8]};
    double teta[3]= {brickData[9],brickData[10],brickData[11]};
    int speedLaw= brickData[12];
    double movingDirection[3]={brickData[13],brickData[14],brickData[15]};
    currentField->info_moving.push_back(speedLaw);
    currentField->info_moving.push_back(movingDirection[0]);
    currentField->info_moving.push_back(movingDirection[1]);
    currentField->info_moving.push_back(movingDirection[2]);
    bool stack=true;
      if(stack == true){
        tempo[0] -= s; tempo[1] -= s; tempo[2] -= s;
    }
    int flag_1 =0;
    int flag_2 =0;
    int flag_3 =0;
   if(tempo[0]==0)
   {
   tempo[0]=s;
   flag_1=1;
   }
   if(tempo[1]==0){
   tempo[1]=s;
   flag_2=1;
   }
   if(tempo[2]==0)
   {
   tempo[2]=s;
   flag_3=1;
   }
    // calculate nb of particles along each direction from target size "s"
    int ni = int(ceil(tempo[0]/s));
    double dx = tempo[0]/ni; ++ni;
    if(flag_1==1){
    ni=ni-1;}
    int nj = int(ceil(tempo[1]/s));
    double dy = tempo[1]/nj; ++nj;
    if(flag_2==1){
    nj=nj-1;}
    int nk = int(ceil(tempo[2]/s));
    double dz = tempo[2]/nk; ++nk;
    if(flag_3==1){
    nk=nk-1;}

 currentField->info_block.push_back((double)c);// tell us if the block is moving, free or fixed
 currentField->info_block.push_back((double)ni);
 currentField->info_block.push_back((double)nj);
 currentField->info_block.push_back((double)nk);
if(c==1 && teta[0]!=0 && teta[1]==0 && teta[2]==0)
{
    currentField->info_block.push_back(teta[0]);
    currentField->info_block.push_back(1);// axe de roation selon x 
}
else if(c==1 && teta[0]==0 && teta[1]!=0 && teta[2]==0)
{
    currentField->info_block.push_back(teta[1]);
    currentField->info_block.push_back(2);// axe de roation selon y 
}
else if(c==1 && teta[0]==0 && teta[1]==0 && teta[2]!=0)
{
    currentField->info_block.push_back(teta[2]);
    currentField->info_block.push_back(3);// axe de roation selon z 
}
else
{
    currentField->info_block.push_back(0);// no roation for the moving 
    currentField->info_block.push_back(0);
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
        }
        break;
        case movingPart :
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
        }
        break;
    }
    return noError;
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
Error readGeometry(std::string filename, Field* currentField, std::vector<double>* volVector)
{
    std::vector<double> posFree, posFixed, posMoving, volVectorFree, volVectorFixed, volVectorMoving;
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
                if(cnt!=N_UL){
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
                    if(readBrick(cube,&inFile, currentField,
                        &posFree, &posMoving, &posFixed,
                        &volVectorFree, &volVectorFixed, &volVectorMoving)==geometryError){
                            return geometryError;
                    }
                }
                else if(buf=="cylin")
                {
                    if(readBrick(cylinder,&inFile, currentField,
                        &posFree, &posMoving, &posFixed,
                        &volVectorFree, &volVectorFixed, &volVectorMoving)==geometryError){
                            return geometryError;
                    }
                }
                else if(buf=="spher")
                {
                    if(readBrick(sphere,&inFile, currentField,
                        &posFree, &posMoving, &posFixed,
                        &volVectorFree, &volVectorFixed, &volVectorMoving)==geometryError){
                            return geometryError;
                    }
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
                    (*volVector)=volVectorFree;
                    currentField->pos=posFree;
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
                            parameter->charactTime=atof(valueArray);
                        if(cnt==9)
                            parameter->c=atof(valueArray);
                        if(cnt==10)
                            parameter->alpha=atof(valueArray);
                        if(cnt==11)
                            parameter->beta=atof(valueArray);
                        if(cnt==12)
                            parameter->epsilon=atof(valueArray);
                        if(cnt==13)
                            parameter->molarMass = atof(valueArray);
                        if (cnt==14)
                            parameter->temperature = atof(valueArray);
                        if (cnt==15)
                            parameter->theta = atof(valueArray);
                        if(cnt==16)
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
                        if(cnt==17)
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
                        if(cnt==18)
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
                        if(cnt==19)
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
                        if(cnt==20)
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
                        if(cnt==21)
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
                        if(cnt==22)
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
                        if(cnt==23)
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
                if(cnt!=N_PARAM){
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
