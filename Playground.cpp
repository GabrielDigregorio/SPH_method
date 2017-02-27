#include "SPH.hpp"

// Structure playground with limited acess to GeneratePlayground function
struct Playground{
    //type of geometry
    std::vector<int> geometry;
    // free geometry
    std::vector<double> geoFreeParam, geoFreeCoord, geoFreeDimen;
    // moving geometry
    std::vector<double> geoMovingParam, geoMovingCoord, geoMovingDimen;
    // fixed geometry
    std::vector<double> geoFixedParam, geoFixedCoord, geoFixedDimen;
};


// Function to fill vectors
void fillVector(std::vector<double> &vect, double A, double B, double C)
{
        vect.push_back(A);
        vect.push_back(B);
        vect.push_back(C);
}


// ReadPlayground: Read the entire .kzr playground and store all data in separate variables
//      Input : filename
//      output: pointer to a structure that contains all data
// !!!SHOULD BE!!! : static Playground *ReadPlayground(const char *filename)
void GeneratePlayground( std::vector<double> &posFree, std::vector<double> &posMoving, std::vector<double> &posFixed, const char *filename)
{
    // open a file geometry.kzr
    std::ifstream infile(filename);
    std::string line;
    std::getline(infile, line);

    // Playground initialisation
    Playground *myPlayground, ptr;
    myPlayground = &ptr;

    // tempory parameters
    int geom;// param[0]= geometry type (cube, cylinder, sphere), 
    double coord[3], dimen[3];
    double param[4];    // param[0]= state of geometry (free, moving, fixed),
                        // param[1]= s (space interval between particles),  
                        // param[2]= % of random position

    // Start Reading The Entire File .kzr
    while (std::getline(infile, line))
    {
        if(line == "#geom") // new geometry environement 
        {
            while (std::getline(infile, line)) // skip white space
            {   
                if(line == "    #brick" || line == "    #cylin" || line == "    #spher")// known geometry
                {
                    if(line == "    #brick")
                        geom = 1; // Cube identifier
                    else if(line == "    #cylin")
                        geom = 2; // Cylinder identifier
                    else if(line == "    #spher")
                        geom = 3; // Sphere identifier
                    
                    myPlayground->geometry.push_back(geom);
                    std::getline(infile, line);

                    // check for parameters
                    if(line == "        #param")
                    {   
                        std::getline(infile, line);
                        param[0] = atof(line.erase(0,10).c_str()); // Status of the geometry
                        std::getline(infile, line);
                        param[1] = atof(line.erase(0,10).c_str()); // Spacing between particles
                        std::getline(infile, line);
                        param[2] = atof(line.erase(0,10).c_str()); // Random particle %
                    }

                    std::getline(infile, line);

                    // check for coordinate
                    if(line == "        #coord")
                    {  
                        std::getline(infile, line);
                        coord[0] = atof(line.erase(0,10).c_str()); // X (center) coordinate of the geometry
                        std::getline(infile, line);
                        coord[1] = atof(line.erase(0,10).c_str()); // Y (center) coordinate of the geometry
                        std::getline(infile, line);
                        coord[2] = atof(line.erase(0,10).c_str()); // Z (center) coordinate of the geometry
                    }

                    std::getline(infile, line); 

                    // check for dimensions
                    if(line == "        #dimen")
                    { 
                        std::getline(infile, line);
                        dimen[0] = atof(line.erase(0,10).c_str()); // Length of the geometry
                        std::getline(infile, line);
                        dimen[1] = atof(line.erase(0,10).c_str()); // Width of the geometry
                        std::getline(infile, line);
                        dimen[2] = atof(line.erase(0,10).c_str()); // Height of the geometry
                    }

                    // Memorise the brick in vectors (separate vector for each status parameter)
                    if(param[0] == 0)
                    {
                        fillVector(myPlayground->geoFreeParam, param[0], param[1], param[2]);
                        fillVector(myPlayground->geoFreeCoord, coord[0], coord[1], coord[2]);
                        fillVector(myPlayground->geoFreeDimen, dimen[0], dimen[1], dimen[2]);
                    }
                    else if(param[0] == 1)
                    {
                        fillVector(myPlayground->geoMovingParam, param[0], param[1], param[2]);
                        fillVector(myPlayground->geoMovingCoord, coord[0], coord[1], coord[2]);
                        fillVector(myPlayground->geoMovingDimen, dimen[0], dimen[1], dimen[2]);
                    }
                    else if(param[0] == 2)
                    {
                        fillVector(myPlayground->geoFixedParam, param[0], param[1], param[2]);
                        fillVector(myPlayground->geoFixedCoord, coord[0], coord[1], coord[2]);
                        fillVector(myPlayground->geoFixedDimen, dimen[0], dimen[1], dimen[2]);
                    }
                    if(1) // put 1 to display value in terminal
                    {
                        std::cout<<"geometry "<<geom<< " , status "<<param[0]<< " , s_spacing "<<param[1]<< " , %random "<<param[2]<<"\n";
                        std::cout<<"coord "<<coord[0]<<" "<<coord[1]<<" "<<coord[2]<<"\n";
                        std::cout<<"dimen "<<dimen[0]<<" "<<dimen[1]<<" "<<dimen[2]<<"\n";
                    }
                }
            }
        }
    }// End Reading The Entire File .kzr



//***********************************************************************
// Generate particles for each geometry !!!SHOULD BE IN AN OTHER FUNCTION!!!
//***********************************************************************

    // For free particles
    for(unsigned int i=0; i<myPlayground->geoFreeCoord.size(); i+=3)
    {
        double o[3] = { myPlayground->geoFreeCoord[i],
                        myPlayground->geoFreeCoord[i+1],
                        myPlayground->geoFreeCoord[i+2]};
        double L[3] = { myPlayground->geoFreeDimen[i],
                        myPlayground->geoFreeDimen[i+1],
                        myPlayground->geoFreeDimen[i+2]};
        double s = myPlayground->geoFreeParam[i+1];
        double r = myPlayground->geoFreeParam[i+2];

        //Generate the geometry for Free particles
        switch (myPlayground->geometry[(i/3)]){

        case 1 : // Cube
            meshcube(o,L,s,posFree, r);
        break;
        case 2 : // Cylinder
            meshcylinder(o,L,s,posFree, r);
        break;
        case 3 : // Sphere
            //meshspherer(o,L,s,posFree, r);
        break;
        }
    }
    // For moving particles
    for(unsigned int i=0; i<myPlayground->geoMovingCoord.size(); i+=3)
    {
        double o[3] = { myPlayground->geoMovingCoord[i],
                        myPlayground->geoMovingCoord[i+1],
                        myPlayground->geoMovingCoord[i+2]};
        double L[3] = { myPlayground->geoMovingDimen[i],
                        myPlayground->geoMovingDimen[i+1],
                        myPlayground->geoMovingDimen[i+2]};
        double s = myPlayground->geoMovingParam[i+1];
        double r = myPlayground->geoMovingParam[i+2];

        //Generate the geometry for Moving particles
        switch (myPlayground->geometry[(i/3)]){

        case 1 : // Cube
            meshcube(o,L,s,posMoving, r);
        break;
        case 2 : // Cylinder
            meshcylinder(o,L,s,posMoving, r);
        break;
        case 3 : // Sphere
            //meshspherer(o,L,s,posMoving, r);
        break;
        }
    }
    // For fixed particles
    for(unsigned int i=0; i<myPlayground->geoFixedCoord.size(); i+=3)
    {
        double o[3] = { myPlayground->geoFixedCoord[i],
                        myPlayground->geoFixedCoord[i+1],
                        myPlayground->geoFixedCoord[i+2]};
        double L[3] = { myPlayground->geoFixedDimen[i],
                        myPlayground->geoFixedDimen[i+1],
                        myPlayground->geoFixedDimen[i+2]};
        double s = myPlayground->geoFixedParam[i+1];
        double r = myPlayground->geoFixedParam[i+2];

        //Generate the geometry for Fixed particles
        switch (myPlayground->geometry[(i/3)]){

        case 1 : // Cube
            meshcube(o,L,s,posFixed, r);
        break;
        case 2 : // Cylinder
            meshcylinder(o,L,s,posFixed, r);
        break;
        case 3 : // Sphere
            //meshspherer(o,L,s,posFixed, r);
        break;
        }
    }

    //return myPlayground;
}

