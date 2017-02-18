#include "sph.h"
#include <cmath>


using namespace std;


/// Naive method to find the neighbors of a given particle at position p(x,y,z)

void neighborAllPair (double p[3], double h, std::vector<double> &pos, std::vector<double> &neighbor)
{
    // memory allocation
    /*pos.reserve(pos.size() + ??? ); GABY : We don't know how many neighbors...*/

    for(int i=0; i<=pos.size(); i=i+3)
    {
        if(pow((p[0]-pos[i]),2)+pow((p[1]-pos[i+1]),2)+pow((p[2]-pos[i+2]),2) <= pow(h,2))
        {
            // the vector neighbor takes the pos vector indices where particles are close to p
            neighbor.push_back(i);
        }

    }

}

/// Advanced method to find the neighbors of a given particle at position p(x,y,z)

// Linked-list algorithm
/*
Inputs:
- pos: array with the particle positions (x,y,z components)
- l: lowest (x,y,z) point in the numerical domain
- u: highest (x,y,z) point in the numerical domain
- kh: support size of the kernel
- values: nonzero values of the incidence matrix (r_ab for the moment...)
- row: row indices of the nonzero values
- column: column indices of the nonzero values
Output:
/
*/
void neighborLinkedList (std::vector<double> &pos,
                         double l[3],
                         double u[3],
                         double kh,
                         std::vector<double> &values,
                         std::vector<int> &row,
                         std::vector<int> &column)
{

    // Box definition
    vector<vector<int> > boxes;
    int nBoxesX = ceil((u[0] - l[0])/kh); // Extra box if non integer quotient
    int nBoxesY = ceil((u[1] - l[1])/kh);
    int nBoxesZ = ceil((u[2] - l[2])/kh);

    int nBoxes = nBoxesX * nBoxesY * nBoxesZ;

    for(int i=0 ; i<nBoxes ; i++)
    {
        vector<int> row;
        boxes.push_back(row);
    }

    //std::cout << nBoxesX << " " << nBoxesY << " " << nBoxesZ << " \n";


    // Sort the particles
    int nPart = pos.size()/3;
    // Box identifier variables
    int boxX;
    int boxY;
    int boxZ;
    double temp;
    for(int i=0 ; i<nPart ; i++)
    {
        temp = (pos[3*i] - l[0])/kh; // Integer division
        boxX = (temp < nBoxesX-1) ? temp : nBoxesX-1;
        temp = (pos[3*i+1] - l[1])/kh;
        boxY = (temp < nBoxesY-1) ? temp : nBoxesY-1;
        temp = (pos[3*i+2] - l[2])/kh;
        boxZ = (temp < nBoxesZ-1) ? temp : nBoxesZ-1;
        // Put the particle identifier in the corresponding box array
        boxes[boxX*nBoxesY*nBoxesZ + boxY*nBoxesZ + boxZ].push_back(i);
        //std::cout << boxX*nBoxesY*nBoxesZ + boxY*nBoxesZ + boxZ << " \n";
    }



    // Search for their neighbors





}

// tree search algorithm
void neighborTree (double p[3], double h, std::vector<double> &pos)
{

}
