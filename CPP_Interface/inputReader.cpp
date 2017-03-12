#include "inputReader.h"

void readBrick(int type, std::ifstream* inFile, Field* currentField, std::vector<int>* sValues){
        std::string buf;
        int cnt=0;
        char valueArray[1024];
        float brickData[N_DATA];
        while(cnt!=N_DATA){
                std::getline(*inFile, buf);
                if(1==sscanf(buf.c_str(),"%*[^=]=%s", valueArray)){
                        brickData[cnt]=atof(valueArray);
                        ++cnt;
                }
                else{continue;}
        }
        int c=(int)brickData[0];
        float s=brickData[1];
        float r=brickData[2];
        double o[3] = {brickData[3],brickData[4],brickData[5]};
        double L[3] = {brickData[6],brickData[7],brickData[8]};
        switch(type){
                case cube :
                        if(c==freePart)
                                meshcube(o, L, s, currentField->posFree, r, true);
                        else if(c==movingPart)
                                meshcube(o, L, s, currentField->posMoving, r, true);
                        else if(c==fixedPart)
                                meshcube(o, L, s, currentField->posFixed, r, true);
                break;
                case cylinder :
                        if(c==freePart)
                                meshcylinder(o, L, s, currentField->posFree, r, true);
                        else if(c==movingPart)
                                meshcylinder(o, L, s, currentField->posMoving, r, true);
                        else if(c==fixedPart)
                                meshcylinder(o, L, s, currentField->posFixed, r, true);
                break;
                case sphere :
                        if(c==freePart)
                                meshsphere(o, L, s, currentField->posFree, r, true);
                        else if(c==movingPart)
                                meshsphere(o, L, s, currentField->posMoving, r, true);
                        else if(c==fixedPart)
                                meshsphere(o, L, s, currentField->posFixed, r, true);
                break;
        }
        sValues->push_back(s); // TO BE USED IF S VALUES ARE NEEDED
}

void readGeometry(std::string filename, Field* currentField){
        std::vector<int> sValues; // TO BE USED IF S VALUES ARE NEEDED
        std::ifstream inFile(filename);
        std::string buf;
        char valueArray[1024];
        int cnt=0;
        std::getline(inFile, buf);
        assert(buf=="#GEOM1" && "Wrong Geometry file type.");
        while(true){ // UNTIL #END_G
                std::getline(inFile, buf);
                buf.erase(std::remove(buf.begin(), buf.end(), ' '),buf.end());
                switch(buf[0]){
                        case '%' :
                                continue;
                        break;
                        case '#' :
                                buf.erase (0,1);
                                if(buf=="domsz"){
                                        while(cnt!=N_UL){
                                                std::getline(inFile, buf);
                                                if(1==sscanf(buf.c_str(),"%*[^=]=%s", valueArray)){
                                                        currentField->u[cnt]=atof(valueArray);
                                                        ++cnt;
                                                }
                                                else
                                                        continue;
                                        }
                                        while(cnt!=N_UL){
                                                std::getline(inFile, buf);
                                                if(1==sscanf(buf.c_str(),"%*[^=]=%s", valueArray)){
                                                        currentField->l[cnt]=atof(valueArray);
                                                        ++cnt;
                                                }
                                                else
                                                        continue;
                                        }
                                }
                                else if(buf=="brick")
                                        readBrick(cube,&inFile, currentField, &sValues);
                                else if(buf=="cylin")
                                        readBrick(cylinder,&inFile, currentField, &sValues);
                                else if(buf=="spher")
                                        readBrick(sphere,&inFile, currentField, &sValues);
                                else if(buf=="END_G"){
                                        return; // REPLACE BY return(0);
                                }
                                else{
                                        std::cout <<"Unknown '"<<buf<<"' identifier."<< "\n";
                                        return; // REPLACE BY return(1);
                                }
                        break;
                        default :
                                continue;
                }
        }
}

void readParameter(std::string filename, Parameter* parameter){
        // open a file geometry.kzr
        std::ifstream inFile(filename);
        std::string buf;
        std::getline(inFile, buf);
        assert(buf=="#PARA1" && "Wrong Parameter file type.");
        while(true){
                std::getline(inFile, buf);
                buf.erase(std::remove(buf.begin(), buf.end(), ' '),buf.end());
                switch(buf[0]){
                        case '%' :
                                continue;
                        break;
                        case '#' :
                                buf.erase (0,1);
                                if(buf=="param"){
                                        int cnt=0;
                                        char valueArray[1024];
                                        while(cnt!=N_PARAM){
                                                std::getline(inFile, buf);
                                                if(1==sscanf(buf.c_str(),"%*[^=]=%s", valueArray)){
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
                                                                parameter->integrationMethod=valueArray;
                                                        if(cnt==9)
                                                                parameter->densityInitMethod=valueArray;
                                                        if(cnt==10)
                                                                parameter->stateEquationMethod=valueArray;
                                                        if(cnt==11)
                                                                parameter->massInitMethod=valueArray;
                                                        if(cnt==12)
                                                                parameter->speedLaw=valueArray;
                                                        ++cnt;
                                                }
                                                else{continue;}
                                        }
                                }
                                else if(buf=="END_F"){
                                        return; // REPLACE BY return(0);
                                }
                                else{
                                        std::cout <<"Unknown '"<<buf<<"' identifier."<< "\n";
                                        return; // REPLACE BY return(1);
                                }
                        break;
                        default :
                                continue;
                }
        }
}
