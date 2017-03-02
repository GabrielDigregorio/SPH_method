#include "SPH.hpp"
#include "Playground.hpp"

/// Constructor : Initialize a 3D matrix and 2D matrix.

Playground :: Playground()
{
    for(int i=0;i<3;++i)
    {
        for(int j=0;j<3;++j)
        {
            DATA.push_back(std::vector< std::vector<double> >());
            DATA[j].push_back(std::vector<double>());
        }
        geometry.push_back(std::vector<int>());
    }
}



/// Function to fill vectors in ReadPlayground function

void Playground :: fillVector(double *A, double *B, double *C, int i)
{
    DATA[i][0].push_back(A[0]); DATA[i][0].push_back(A[1]); DATA[i][0].push_back(A[2]);
    DATA[i][1].push_back(B[0]); DATA[i][1].push_back(B[1]); DATA[i][1].push_back(B[2]);
    DATA[i][2].push_back(C[0]); DATA[i][2].push_back(C[1]); DATA[i][2].push_back(C[2]);
}



/// ReadPlayground: Read the entire .kzr playground and store all data 
    //                 in separate variables within a structure
    //                 (NO GENERATION OF PARTICLES HERE)
    //      Input : filename
    //      output: structure that contains all data

void Playground :: ReadPlayground(const char *filename)
{
    // open a file geometry.kzr
    std::ifstream infile(filename);
    std::string line;
    std::getline(infile, line);

    // tempory parameters
    int geom;// geometry type (cube, cylinder, sphere), 
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
                    if     (line == "    #brick")
                        geom = 1; // Cube identifier
                    else if(line == "    #cylin")
                        geom = 2; // Cylinder identifier
                    else if(line == "    #spher")
                        geom = 3; // Sphere identifier
                    
                    std::getline(infile, line);

                    // check for parameters
                    if(line == "        #param")
                    {   
                        std::getline(infile, line);
                        param[0] = atof(line.erase(0,10).c_str()); // Status of the geometry

                        //if(param[0] == 0)
                            //geoFree.push_back(geom); // Free identifier
                        //else if(param[0] == 1)
                            //geoMoving.push_back(geom); // Moving identifier
                        //else if(param[0] == 2)
                            //geoFixed.push_back(geom); // Fixed identifier

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
                    fillVector(param, coord, dimen, param[0]);
                    geometry[param[0]].push_back(geom);

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

}



/// GeneratePlayground: Generate all particles in all geometries from structure Playground
    //      Input : posFree, posMoving, posFixed, filename
    //      output: filled posFree, posMoving, posFixed by structure Playground

void Playground :: GeneratePlayground( std::vector<double> &posFree, std::vector<double> &posMoving, std::vector<double> &posFixed)
{
    //Stack all geometries
    bool stack = true;

    for(int c=0; c<3; ++c)
    {
        // For free particles
        for(int i=0; i<DATA[c][1].size(); i+=3)
        {
            double o[3] = { DATA[c][1][i],
                            DATA[c][1][i+1],
                            DATA[c][1][i+2]};
            double L[3] = { DATA[c][2][i],
                            DATA[c][2][i+1],
                            DATA[c][2][i+2]};
            double s = DATA[c][0][i+1];
            double r = DATA[c][0][i+2];

            //Generate the geometry for Free particles
            switch (geometry[c][(i/3)]){
            case 1 : // Cube
                meshcube(o, L, s, posFree, r, stack);
            break;
            case 2 : // Cylinder
                meshcylinder(o, L, s, posFree, r, stack);
            break;
            case 3 : // Sphere
                meshsphere(o, L, s, posFree, r, stack);
            break;
            }
        }
    }

}