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

void Playground :: fillVector(std::vector<double> data, int i)
{
    DATA[i][0].push_back(data[0]); DATA[i][0].push_back(data[1]); DATA[i][0].push_back(data[2]);
    DATA[i][1].push_back(data[3]); DATA[i][1].push_back(data[4]); DATA[i][1].push_back(data[5]);
    DATA[i][2].push_back(data[6]); DATA[i][2].push_back(data[7]); DATA[i][2].push_back(data[8]);
}



/// ReadPlayground: Read the entire .kzr playground and store all data 
    //                 in separate variables within a structure
    //                 (NO GENERATION OF PARTICLES HERE)
    //      Input : filename .kzr
    //      output: structure that contains all data

void Playground :: ReadPlayground(const char *filename)
{
    // open a file geometry.kzr
    std::ifstream infile(filename);
    std::string line;
    std::getline(infile, line);

    int geom;

    // Start Reading The Entire File .kzr
    while (std::getline(infile, line))
    {
        if(line == "#geom") // new geometry environement 
        {
            while (std::getline(infile, line)) // skip white space
            {   
                // known geometry 
                if(line == "    #brick" || line == "    #cylin" || line == "    #spher")
                {
                    if     (line == "    #brick")
                        geom = 1; // Cube identifier
                    else if(line == "    #cylin")
                        geom = 2; // Cylinder identifier
                    else if(line == "    #spher")
                        geom = 3; // Sphere identifier
                    
                    std::getline(infile, line);

                    // Load data for a given geometry
                    for(int i=0; i<nbr_data; ++i)
                    {
                        std::getline(infile, line);
                        data[i]=atof(line.erase(0,10).c_str());
                        if((i+1)%3 == 0)
                            std::getline(infile, line);
                    }

                    // Memorise the brick in vectors (separate vector for each status parameter)
                    fillVector(data, data[0]);
                    geometry[data[0]].push_back(geom);

                    if(1) // put 1 to display value in terminal
                    {
                        std::cout<<"geometry "<<geom<< " , status "<<data[0]<<
                                   " , s_spacing "<<data[1]<< " , %random "<<data[2]<<"\n";
                        std::cout<<"coord "<<data[3]<<" "<<data[4]<<" "<<data[5]<<"\n";
                        std::cout<<"dimen "<<data[6]<<" "<<data[7]<<" "<<data[8]<<"\n";
                    }
                }
            }
        }
    }// End Reading The Entire File .kzr

}



/// GeneratePlayground: Generate all particles in all geometries from structure Playground
    //      Input : posFree, posMoving, posFixed, filename
    //      output: filled posFree, posMoving, posFixed by structure Playground

void Playground :: GeneratePlayground(  std::vector<double> &posFree, 
                                        std::vector<double> &posMoving, 
                                        std::vector<double> &posFixed)
{
    //Stack all geometries
    bool stack = true;

    for(int c=0; c<3; ++c)
    {
        // For free particles
        for(int i=0; i<DATA[c][1].size(); i+=3)
        {
            double o[3] = { DATA[c][1][i], DATA[c][1][i+1], DATA[c][1][i+2]};
            double L[3] = { DATA[c][2][i], DATA[c][2][i+1], DATA[c][2][i+2]};
            double s    =   DATA[c][0][i+1];
            double r    =   DATA[c][0][i+2];

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