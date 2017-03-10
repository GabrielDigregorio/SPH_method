#include "Main.h"
#include "Physics.h"

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
                         std::vector<std::vector<int> > &neighborsAll,
                         std::vector<std::vector<double> > &kernelGradientsAll)
{
    double kh2 = kh*kh;
    // For each particle, browse all other particles and compute the distance
    for(unsigned int i=0; i<pos.size()/3; i++){
        for(unsigned int j=0; j<pos.size()/3; j++){
            double r2 = distance(pos, i, j);
            if( r2 < kh2 ){
                neighborsAll[i].push_back(j); // The neighbor ID
                kernelGradientsAll[i].push_back(2.0); // The kernel grandient with this neighbor
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
void neighborLinkedList(std::vector<double> &pos,
                         double l[3],
                         double u[3],
                         double kh,
                         std::vector<std::vector<int> > &neighborsAll,
                         std::vector<std::vector<double> > &kernelGradientsAll)
{
    double kh2 = kh*kh; // More efficient to compare distance^2

    // Box definition
    std::vector<std::vector<int> > boxes;
    int nBoxesX = ceil((u[0] - l[0])/kh); // Extra box if non integer quotient
    int nBoxesY = ceil((u[1] - l[1])/kh);
    int nBoxesZ = ceil((u[2] - l[2])/kh);

    int nBoxes = nBoxesX * nBoxesY * nBoxesZ;

    for(int i=0 ; i<nBoxes ; i++){
        std::vector<int> boxContent;
        boxes.push_back(boxContent);
    }

    // Sort the particles
    int nPart = pos.size()/3;
    // Box identifier variables
    int boxX;
    int boxY;
    int boxZ;
    double temp;
    for(int i=0 ; i<nPart ; i++){
        temp = (pos[3*i] - l[0])/kh; // Integer division
        boxX = (temp < nBoxesX-1) ? temp : nBoxesX-1;
        temp = (pos[3*i+1] - l[1])/kh;
        boxY = (temp < nBoxesY-1) ? temp : nBoxesY-1;
        temp = (pos[3*i+2] - l[2])/kh;
        boxZ = (temp < nBoxesZ-1) ? temp : nBoxesZ-1;
        // Put the particle identifier in the corresponding box array
        boxes[boxX + boxY*nBoxesX + boxZ*nBoxesX*nBoxesY].push_back(i);
    }

    // Search for their neighbors
    int particleID;
    // Spans the boxes
    for(int box=0 ; box<nBoxes ; box++){
        // Determines the list of surronding boxes (boundaries -> not trivial)
        std::vector<int> surrBoxes;
        surroundingBoxes(box, nBoxesX, nBoxesY, nBoxesZ, surrBoxes);
        // Spans the particles in the box
        for(unsigned int part=0 ; part<boxes[box].size() ; part++){
            particleID = boxes[box][part];
            // Spans the surrounding boxes
            for(unsigned int surrBox = 0 ; surrBox < surrBoxes.size() ; surrBox++){
                // Spans the higher index particles in the box (symmetry)
                for(unsigned int i=0 ; i<boxes[surrBoxes[surrBox]].size() ; i++){
                    int potNeighborID = boxes[surrBoxes[surrBox]][i];
                    double r2 = distance(pos, particleID, potNeighborID);
                    if(r2<kh2){
                        neighborsAll[particleID].push_back(potNeighborID);
                        kernelGradientsAll[particleID].push_back(2.0);
                    }
                }
            }
        }
    }
}


/* ----------------------- UPDATE ----------------------*/

// Creates a mesh to sort the particles and gives the box adjacent relations
void boxMesh(double l[3], double u[3], double kh,
            std::vector<std::vector<int> > &boxes,
            std::vector<std::vector<int> > &surrBoxesAll){
  // Determination of the number of boxes in each direction
  int nBoxesX = ceil((u[0] - l[0])/kh); // Extra box if non integer quotient
  int nBoxesY = ceil((u[1] - l[1])/kh);
  int nBoxesZ = ceil((u[2] - l[2])/kh);
  int nBoxes = nBoxesX * nBoxesY * nBoxesZ;

  // Determines the neighboring relations (DOES NOT CREATE THE BOXES)
  for(int box=0 ; box<nBoxes ; box++){
      std::vector<int> boxContent;
      std::vector<int> surrBoxes;
      surroundingBoxes(box, nBoxesX, nBoxesY, nBoxesZ, surrBoxes); // Fills the list
      boxes.push_back(boxContent); // Add the (empty) box vector
      surrBoxesAll.push_back(surrBoxes); // Add the (filled) surrounding box list.
  }
  return;
}

// Sorts the particles into cubic boxes
void sortParticles(std::vector<double> &pos, double l[3], double u[3], double kh,
                   std::vector<std::vector<int> > &boxes){
    // Determination of the number of boxes in each direction
    int nBoxesX = ceil((u[0] - l[0])/kh); // Extra box if non integer quotient
    int nBoxesY = ceil((u[1] - l[1])/kh);
    int nBoxesZ = ceil((u[2] - l[2])/kh);
    // Box identifier variables
    int boxX; int boxY; int boxZ;
    double temp;
    for(int i=0 ; i<pos.size()/3 ; i++){
        temp = (pos[3*i] - l[0])/kh; // Integer division
        boxX = (temp < nBoxesX-1) ? temp : nBoxesX-1;
        temp = (pos[3*i+1] - l[1])/kh;
        boxY = (temp < nBoxesY-1) ? temp : nBoxesY-1;
        temp = (pos[3*i+2] - l[2])/kh;
        boxZ = (temp < nBoxesZ-1) ? temp : nBoxesZ-1;
        // Put the particle identifier in the corresponding box array
        boxes[boxX + boxY*nBoxesX + boxZ*nBoxesX*nBoxesY].push_back(i);
    }
}

/* Searches the neighbors of a given particle in the surrounding boxes
Fills the neighbors/kernelGradients vectors with the neighbors and the associated
values of the kernel gradient for the given particleID.
*/
void findNeighbors(int particleID, std::vector<double> &pos, double kh2,
                   std::vector<std::vector<int> > &boxes,
                   std::vector<int> &surrBoxes,
                   std::vector<int> &neighbors,
                   std::vector<double> &kernelGradients){
    // Spans the surrounding boxes
    for(unsigned int surrBox = 0 ; surrBox < surrBoxes.size() ; surrBox++){
        // Spans the particles in the box (all particles!)
        for(unsigned int i=0 ; i<boxes[surrBoxes[surrBox]].size() ; i++){
            int potNeighborID = boxes[surrBoxes[surrBox]][i];
            double r2 = distance(pos, particleID, potNeighborID);
            if(r2<kh2){
                neighbors.push_back(potNeighborID);
                double kernelGradient = 1.0; // TO CHANGE !!!!!!!
                kernelGradients.push_back(kernelGradient);
            }
        }
    }
}

// Clear the boxes content
void boxClear(std::vector<std::vector<int> > &boxes){
    for(int i=0 ; i<boxes.size() ; i++)
        boxes[i].clear();
}


// Gives the list of the surrounding boxes
void surroundingBoxes(int box, int nBoxesX, int nBoxesY, int nBoxesZ, std::vector<int> &surrBoxes){
    int index_x, index_y, index_z;
    index_z = box/(nBoxesX*nBoxesY);
    index_y = (box-index_z*nBoxesX*nBoxesY)/nBoxesX;
    index_x = box-index_z*nBoxesX*nBoxesY-index_y*nBoxesX;

    std::vector<int> tab(6, 1); // Initializes to 1
    std::vector<int> value(9, 0); // Initialized to 0
    value[0]=-1; value[2]=1;
    value[3]=-nBoxesX; value[5]=nBoxesX;
    value[6]=-nBoxesX*nBoxesY; value[8]=nBoxesX*nBoxesY;

    // Filling of the tab vector.
    if (index_x>0){tab[0]=0;}
    if (index_x<nBoxesX-1){tab[1]=2;}
    if (index_y>0){tab[2]=0;}
    if (index_y<nBoxesY-1){tab[3]=2;}
    if (index_z>0){tab[4]=0;}
    if (index_z<nBoxesZ-1){tab[5]=2;}

    // Finding the neighbors
    for (int k = tab[4]; k <= tab[5]; k++){
        for (int j = tab[2]; j <= tab[3]; j++){
            for (int i = tab[0]; i <= tab[1]; i++){
                surrBoxes.push_back(box+value[i]+value[j+3]+value[k+6]);
            }
        }
    }
    return;
}


// Gives the distance between two particles
double distance(std::vector<double> &pos, int partA, int partB){
    return (pos[partA*3]-pos[partB*3])*(pos[partA*3]-pos[partB*3])
             + (pos[partA*3+1]-pos[partB*3+1])*(pos[partA*3+1]-pos[partB*3+1])
             + (pos[partA*3+2]-pos[partB*3+2])*(pos[partA*3+2]-pos[partB*3+2]);

}
