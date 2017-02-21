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
    int particleID;
    // Spans the boxes
    for(int box=0 ; box<nBoxes ; box++)
    {
        // Determines the list of surronding boxes (boundaries -> not trivial)
        std::vector<int> surrBoxes;
        surroundingBoxes(box, nBoxesX, nBoxesY, nBoxesZ, surrBoxes);
        // Spans the particles in the box
        for(int part=0 ; part<boxes[box].size() ; part++)
        {
            particleID = boxes[box][part];
            // Spans the surrounding boxes
            for(int surrBox = 0 ; surrBox < surrBoxes.size() ; surrBox++)
            {
                // Spans the higher index particles in the box (symmetry)
                for(int i=0 ; i<boxes[surrBox].size() ; i++)
                {
                    int potNeighborID = boxes[surrBox][i];
                    if(potNeighborID >= particleID)
                    {
                        double r = distance(pos, particleID, potNeighborID);
                        if(r<kh)
                        {
                            std::cout << r << "\t: ";
                            std::cout << particleID << " " << potNeighborID << "\n";
                            values.push_back(r); // The distance bewteen the two found neighbors
                            row.push_back(particleID); // The one we search the neighbors of
                            column.push_back(potNeighborID); // The neighbor we have just found
                        }
                    }
                }
            }
        }
    }

}

// Gives the list of the surrounding boxes
void surroundingBoxes(int box, int nBoxesX, int nBoxesY, int nBoxesZ, std::vector<int> &surrBoxes)
{
    surrBoxes.push_back(box);
    // TO COMPLETE !!!
    /* NAIVE WAY !
    // in current plane
    north_neighbors = box - nBoxesX
    northeast_neighbors = box - nBoxesX + 1
    east_neighbors = box + 1
    southeast_neighbors = box + nBoxesX + 1
    south_neighbors = box + nBoxesX
    southwest_neighbors = box + nBoxesX - 1
    west_neighbors = box - 1
    northwest_neighbors = box - nBoxesX - 1

    // in plane z = -1
    front_neighbors = - nBoxesY*nBoxesX
    front_north_neighbors = box - nBoxesX - nBoxesY*nBoxesX
    front_northeast_neighbors = box - nBoxesX + 1 - nBoxesY*nBoxesX
    front_east_neighbors = box + 1 - nBoxesY*nBoxesX
    front_southeast_neighbors = box + nBoxesX + 1 - nBoxesY*nBoxesX
    front_south_neighbors = box + nBoxesX - nBoxesY*nBoxesX
    front_southwest_neighbors = box + nBoxesX - 1 - nBoxesY*nBoxesX
    front_west_neighbors = box - 1 - nBoxesY*nBoxesX
    front_northwest_neighbors = box - nBoxesX - 1 - nBoxesY*nBoxesX

    // in plane z = 1
    back_neighbors = + nBoxesY*nBoxesX
    back_north_neighbors = box - nBoxesX + nBoxesY*nBoxesX
    back_northeast_neighbors = box - nBoxesX + 1 + nBoxesY*nBoxesX
    back_east_neighbors = box + 1 + nBoxesY*nBoxesX
    back_southeast_neighbors = box + nBoxesX + 1 + nBoxesY*nBoxesX
    back_south_neighbors = box + nBoxesX + nBoxesY*nBoxesX
    back_southwest_neighbors = box + nBoxesX - 1 + nBoxesY*nBoxesX
    back_west_neighbors = box - 1 + nBoxesY*nBoxesX
    back_northwest_neighbors = box - nBoxesX - 1 + nBoxesY*nBoxesX

    // Check all cases
     if(z_bondaries = front)
     {
        neighbors in current plane
        neighbors in plane z = +1

        if(y_bondaries = top)
        {
            if(x_bondaries = left)
            {

            }
            else if(x_bondaries = right)
            {

            }
            else
            {

            }
        }
        else if (y_bondaries = bottom)
        {
            if(x_bondaries = left)
            {

            }
            else if(x_bondaries = right)
            {

            }
            else
            {

            }
        }
        else
        {

        }

     }
     else if (z_bondaries = back)
     {
        neighbors in plane z = -1
        neighbors in current plane
        
        if(y_bondaries = top)
        {
            if(x_bondaries = left)
            {

            }
            else if(x_bondaries = right)
            {
                
            }
            else
            {

            }
        }
        else if (y_bondaries = bottom)
        {
            if(x_bondaries = left)
            {

            }
            else if(x_bondaries = right)
            {

            }
            else
            {

            }
        }
        else
        {

        }

     }
     else
     {
        neighbors in plane z = -1
        neighbors in current plane
        neighbors in plane z = +1
        
        if(y_bondaries = top)
        {
            if(x_bondaries = left)
            {

            }
            else if(x_bondaries = right)
            {
                
            }
            else
            {

            }
        }
        else if (y_bondaries = bottom)
        {
            if(x_bondaries = left)
            {

            }
            else if(x_bondaries = right)
            {

            }
            else
            {

            }

        }
        else
        {

        }

     }

    */
    return;
}

// Gives the distance between two particles
double distance(std::vector<double> pos, int partA, int partB)
{
    return sqrt( pow(pos[partA*3]-pos[partB*3],2)
               + pow(pos[partA*3+1]-pos[partB*3+1],2)
               + pow(pos[partA*3+2]-pos[partB*3+2],2));
}




// tree search algorithm
void neighborTree (double p[3], double h, std::vector<double> &pos)
{

}
