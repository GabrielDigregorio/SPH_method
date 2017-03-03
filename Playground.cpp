#include "./Headers/SPH.hpp"
#include "./Headers/Playground.hpp"

/// Private Constructor : Initialize a 3D matrix and 2D matrix.

Playground :: Playground( bool s)
{
    screen = s;
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



/// Private Function to fill vectors in ReadPlayground function

void Playground :: fillVector(std::vector<double> data, int i)
{
    DATA[i][0].push_back(data[0]); DATA[i][0].push_back(data[1]); DATA[i][0].push_back(data[2]);
    DATA[i][1].push_back(data[3]); DATA[i][1].push_back(data[4]); DATA[i][1].push_back(data[5]);
    DATA[i][2].push_back(data[6]); DATA[i][2].push_back(data[7]); DATA[i][2].push_back(data[8]);
}



/// Public ReadPlayground: Read the entire .kzr playground and store all data 
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

    // Start Reading FLUID in File .kzr
    while (line != "#FLUID")
        std::getline(infile, line);
    
    // Read Parameters of the fluid
    for(int i=0; i<6; ++i)
    {
        std::getline(infile, line);
        param.push_back( atof(line.erase(0,8).c_str()) );
    }

    // Read method of the solver (Euler or RungeKutta)
    std::getline(infile, line);
    method = line.erase(0,8).c_str();

    if(screen)
    {
        std::cout<<"\n" << "Fluid: " <<", cst1="<<param[0]<<", cst2="<<param[1]<<", cst3="<<param[2]
                    << ", cst4=" <<param[3]<<", cst5="<<param[4]<<", cst6="<<param[5] <<"\n";
        std::cout<<"Method: "<< method <<"\n"; 
    }
          
    int geom;

    // Start Reading GEOM in File .kzr
    while (line != "#GEOM")
        std::getline(infile, line);

    for(int i=0; i<3; ++i)
    {
        std::getline(infile, line);
        l[i] = atof(line.erase(0,3).c_str());
        std::getline(infile, line);
        u[i] = atof(line.erase(0,3).c_str());
    }

    if(screen)
        std::cout<<"\n"<< "domain: " <<"lx="<<l[0]<<", ly="<<l[1]<<", lz="<<l[2] << 
                            " and "<< " ux=" <<u[0]<<", uy="<<u[1]<<", uz="<<u[2] <<"\n\n";
    
    // For each geometry
    while (std::getline(infile, line)) // check for the end of the file
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

            if(screen) // put 1 to display value in terminal
            { 
                std::cout<<"geometry "<<geom<< " , status "<<data[0]<<
                            " , s_spacing "<<data[1]<< " , %random "<<data[2]<<"\n";
                std::cout<<"coord "<<data[3]<<" "<<data[4]<<" "<<data[5]<<"\n";
                std::cout<<"dimen "<<data[6]<<" "<<data[7]<<" "<<data[8]<<"\n \n";
            }
        }
    }// End Reading The Entire File .kzr

}



/// Public GeneratePlayground: Generate all particles in all geometries from structure Playground
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

// Return the fluid parameters
std::vector<double> Playground :: GetParam()
{
    return param;
}


// Return the method for the solver
std::string Playground :: GetMethod()
{
    return method;
}


// Return the lower coordinate or the upper coordinate of the domain
std::vector<double> Playground :: GetDomain(bool dom)
{
    if(dom == false)
        return l;
    else if (dom == true)
        return u;
}
