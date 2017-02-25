#include "SPH.h"
#include <cmath>

/// Naive method to find the neighbors of all particles.
// All-pair search algorithm
/*
Inputs:
- pos: array with the particle positions (x,y,z components)
- kh: support size of the kernel
- values: nonzero values of the incidence matrix (r_ab for the moment...)
- row: row indices of the nonzero values
- column: column indices of the nonzero values
Output:
/
*/
void neighborAllPair (std::vector<double> &pos,
                         double kh,
                         std::vector<double> &values,
                         std::vector<int> &row,
                         std::vector<int> &column)
{
    double kh2 = pow(kh,2);

    // For each particle, browse all other particles and compute the distance
    for(int i=0; i<pos.size(); i=i+3)
    {
        for(int j=i; j<pos.size(); j=j+3)
        {
            double r2 = distance(pos, i/3, j/3);
            if( r2 < kh2 )
            {
                //std::cout << r2 << "\t: ";
                //std::cout << i << " " << j << "\n";
                values.push_back(r2); // The distance bewteen the two found neighbors
                row.push_back(i); // The one we search the neighbors of
                column.push_back(j); // The neighbor we have just found
            }
        }
    }

}



/// Advanced method to find the neighbors of all particles.

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
    double kh2 = pow(kh,2); // More efficient to compare distance^2

    // Box definition
    std::vector<std::vector<int> > boxes;
    int nBoxesX = ceil((u[0] - l[0])/kh); // Extra box if non integer quotient
    int nBoxesY = ceil((u[1] - l[1])/kh);
    int nBoxesZ = ceil((u[2] - l[2])/kh);

    int nBoxes = nBoxesX * nBoxesY * nBoxesZ;

    for(int i=0 ; i<nBoxes ; i++)
    {
        std::vector<int> boxContent;
        boxes.push_back(boxContent);
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
                        double r2 = distance(pos, particleID, potNeighborID);
                        if(r2<kh2)
                        {
                            //std::cout << r2 << "\t: ";
                            //std::cout << particleID << " " << potNeighborID << "\n";
                            values.push_back(r2); // The distance bewteen the two found neighbors
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
    int nx=nBoxesX;
    int ny=nBoxesY;
    int nz=nBoxesZ;
    int index_x;
    int index_y;
    int index_z;
    index_z=box/(nBoxesX*nBoxesY);
    index_y=(box-index_z*nBoxesX*nBoxesY)/nBoxesX;
    index_x=box-index_z*nBoxesX*nBoxesY-index_y*nBoxesX;
    std::vector<int> tab(6, 1);
    std::vector<int> value(9, 0);
    value[0]=-1;value[1]=0;value[2]=1;value[3]=-nx;value[4]=0;value[5]=nx;value[6]=-nx*ny;value[7]=0;value[8]=nx*ny;
    // remplissage du vecteur tab.
    if (index_x>0)
    {
      tab[0]=0;
    }
    if (index_x<nx-1)
    {
      tab[1]=2;
    }
    if (index_y>0)
    {
      tab[2]=0;
    }
    if (index_y<ny-1)
    {
      tab[3]=2;
    }
    if (index_z>0)
    {
      tab[4]=0;
    }
    if (index_y<ny-1)
    {
      tab[5]=2;
    }

    for (int i = tab[0]; i < tab[1]; i++)
    {

        for (int j = tab[2]; j < tab[3]; j++)
        {

            for (int k = tab[4]; k < tab[5]; k++)
            {

                        // given i j k , we have the block to push
                        // i: 0->-1
                        //    1->0
                        //    2->+1
                        // j: 0->-nx
                        //    1->0
                        //    2->+nx
                        // k: 0->-nx*ny
                        //    1->0
                        //    2->+nx*ny
                        // surrBoxes.push_back(box+...);
                        surrBoxes.push_back(box+value[i]+value[j+3]+value[k+6]);

            }
        }
    }
    return;
}


// Gives the distance between two particles
double distance(std::vector<double> pos, int partA, int partB)
{
    return pow(pos[partA*3]-pos[partB*3],2)
               + pow(pos[partA*3+1]-pos[partB*3+1],2)
               + pow(pos[partA*3+2]-pos[partB*3+2],2);
}


// tree search algorithm
void neighborTree (double p[3], double h, std::vector<double> &pos)
{

}
