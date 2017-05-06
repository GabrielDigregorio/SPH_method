#include "Main.h"
#include "Physics.h"
#include "Tools.h"
#include "Structures.h"

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
void neighborAllPair (std::vector<double> (&pos)[3],
                        double kh,
                        std::vector<std::vector<int> > &neighborsAll,
                        std::vector<std::vector<double> > &kernelGradientsAll,
                        Kernel myKernel)
{
    double kh2 = kh*kh;
    double r;
    double r2;
    double direction;
    double currentKernelGradientMag;
    // For each particle, browse all other particles and compute the distance
    for(unsigned int i=0; i<pos[0].size(); i++){
        for(unsigned int j=0; j<pos[0].size(); j++){
            r2 = distance(pos, i, j);
            if( r2 < kh2 && i != j){
                // Neighbor saving
                neighborsAll[i].push_back(j);
                // Kernel gradient saving
                r = sqrt(r2);
                currentKernelGradientMag = gradWab(r, kh, myKernel);
                for(int coord=0 ; coord<3 ; coord++){
                    direction = (pos[coord][i]-pos[coord][j]) / r;
                    kernelGradientsAll[i].push_back(direction * currentKernelGradientMag);
                }
            }
        }
    }
}



/// Advanced method to find the neighbors of all particles.

// Creates a mesh to sort the particles and gives the box adjacent relations
void boxMesh(double l[3], double u[3], double boxSize,
            std::vector<std::vector<int> > &boxes,
            std::vector<std::vector<int> > &surrBoxesAll){
  // Determination of the number of boxes in each direction
  int nBoxesX = ceil((u[0] - l[0])/boxSize); // Extra box if non integer quotient
  int nBoxesY = ceil((u[1] - l[1])/boxSize);
  int nBoxesZ = ceil((u[2] - l[2])/boxSize);
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


// Determines the box size depending on the integrationMethod
double boxSizeCalc(double kh, IntegrationMethod method){
    double multiplicator = (method)? 1.1 : 1.0;
    return kh * multiplicator;
}

// Sorts the particles into cubic boxes
void sortParticles(std::vector<double> (&pos)[3], double l[3], double u[3], double boxSize,
                   std::vector<std::vector<int> > &boxes){
    // Start from scratch
    boxClear(boxes);
    // Determination of the number of boxes in each direction
    int nBoxesX = ceil((u[0] - l[0])/boxSize); // Extra box if non integer quotient
    int nBoxesY = ceil((u[1] - l[1])/boxSize);
    int nBoxesZ = ceil((u[2] - l[2])/boxSize);
    // Box identifier variables
    int boxX; int boxY; int boxZ;
    double temp;
    // Note: no OpenMP because the push_back would require "critical" section -> not efficient at all (has been tested)
    //      but not a big problem, this loop is not too much time consuming
    for(int i=0 ; i<pos[0].size() ; i++){
        // Box coordinate along X
        temp = (pos[0][i] - l[0])/boxSize; // Integer division
        if(temp < 0){boxX = 0;}
        else{boxX = (temp < nBoxesX-1) ? temp : nBoxesX-1;}
        // Box coordinate along Y
        temp = (pos[1][i] - l[1])/boxSize;
        if(temp < 0){boxY = 0;}
        else{boxY = (temp < nBoxesY-1) ? temp : nBoxesY-1;}
        // Box coordinate along Z
        temp = (pos[2][i] - l[2])/boxSize;
        if(temp < 0){boxZ = 0;}
        else{boxZ = (temp < nBoxesZ-1) ? temp : nBoxesZ-1;}
        // Put the particle identifier in the corresponding box array
        boxes[boxZ + boxY*nBoxesZ + boxX*nBoxesZ*nBoxesY].push_back(i);
    }
}

// Overload with "optimization" -> useless with vectors, maybe useful with lists...
void sortParticles(std::vector<double> (&pos)[3], double l[3], double u[3], double boxSize,
                   std::vector<std::vector<int> > &boxes, bool toOptimize){

    if(toOptimize){
        // Determination of the number of boxes in each direction
        int nBoxesX = ceil((u[0] - l[0])/boxSize); // Extra box if non integer quotient
        int nBoxesY = ceil((u[1] - l[1])/boxSize);
        int nBoxesZ = ceil((u[2] - l[2])/boxSize);

        int boxNumber = 0;
        for(int currBoxX=0 ; currBoxX<nBoxesX ; currBoxX++){
            for(int currBoxY=0 ; currBoxY<nBoxesY ; currBoxY++){
                for(int currBoxZ=0 ; currBoxZ<nBoxesZ ; currBoxZ++){
                    int size = boxes[boxNumber].size();
                    for(int i=0 ; i<size ;){
                        int particleID = boxes[boxNumber][i];
                        int boxX, boxY, boxZ;
                        // Box pre-coordinates along X
                        boxX = (pos[0][particleID] - l[0])/boxSize; // Integer division
                        boxY = (pos[1][particleID] - l[1])/boxSize;
                        boxZ = (pos[2][particleID] - l[2])/boxSize;
                        if(boxX==currBoxX && boxY==currBoxY && boxZ==currBoxZ){i++;}
                        else{
                            // Box coordinate along X
                            if(boxX < 0){boxX = 0;}
                            else{boxX = (boxX < nBoxesX-1) ? boxX : nBoxesX-1;}
                            // Box coordinate along Y
                            if(boxY < 0){boxY = 0;}
                            else{boxY = (boxY < nBoxesY-1) ? boxY : nBoxesY-1;}
                            // Box coordinate along Z
                            if(boxZ < 0){boxZ = 0;}
                            else{boxZ = (boxZ < nBoxesZ-1) ? boxZ : nBoxesZ-1;}
                            // Put the particle identifier in the corresponding box array
                            boxes[boxNumber].erase(boxes[boxNumber].begin()+i);
                            size--;
                            boxes[boxZ + boxY*nBoxesZ + boxX*nBoxesZ*nBoxesY].push_back(particleID);
                        }
                    }
                    boxNumber++;
                }
            }
        }
    }
    else{sortParticles(pos, l, u, boxSize, boxes);}
}



/* Searches the neighbors of a given particle in the surrounding boxes
Fills the neighbors/kernelGradients vectors with the neighbors and the associated
values of the kernel gradient for the given particleID.
*/
void findNeighbors(int particleID, std::vector<double> (&pos)[3], double kh,
                    std::vector<std::vector<int> > &boxes,
                    std::vector<int> &surrBoxes,
                    std::vector<int> &neighbors,
                    std::vector<double> &kernelGradients,
                    Kernel myKernel){
    double kh2 = kh*kh;
    double r;
    double currentKernelGradientMag;
    double direction;
    // Spans the surrounding boxes
    for(unsigned int surrBox = 0 ; surrBox < surrBoxes.size() ; surrBox++){
        // Spans the particles in the box (all particles!)
        for(unsigned int i=0 ; i<boxes[surrBoxes[surrBox]].size() ; i++){
            int potNeighborID = boxes[surrBoxes[surrBox]][i];
            double r2 = distance(pos, particleID, potNeighborID);
            if(r2<kh2 && particleID != potNeighborID){
                // Neighbor saving
                neighbors.push_back(potNeighborID);
                // Kernel gradient saving
                r = sqrt(r2);
                currentKernelGradientMag = gradWab(r, kh, myKernel);
                for(int coord=0 ; coord<3 ; coord++){
                    direction = (pos[coord][particleID]-pos[coord][potNeighborID]) / r;
                    kernelGradients.push_back(direction * currentKernelGradientMag);
                }
            }
        }
    }
}

/* Overload with tabulated values : USELESS !!*/
void findNeighbors(int particleID, std::vector<double> (&pos)[3], double kh,
                    std::vector<std::vector<int> > &boxes,
                    std::vector<int> &surrBoxes,
                    std::vector<int> &neighbors,
                    std::vector<double> &kernelGradients,
                    Kernel myKernel,
                    std::vector<double> &kernelGradientsSamples,
                    int resolution){
    double kh2 = kh*kh;
    double r;
    double currentKernelGradientMag;
    double direction; // will contain (x_a-x_b)/r_ab (or y, z)
    // Spans the surrounding boxes
    for(unsigned int surrBox = 0 ; surrBox < surrBoxes.size() ; surrBox++){
        // Spans the particles in the box (all particles!)
        for(unsigned int i=0 ; i<boxes[surrBoxes[surrBox]].size() ; i++){
            int potNeighborID = boxes[surrBoxes[surrBox]][i];
            double r2 = distance(pos, particleID, potNeighborID);
            if(r2<kh2 && particleID != potNeighborID){
                // Neighbor saving
                neighbors.push_back(potNeighborID);
                // Kernel gradient saving
                r = sqrt(r2);
                currentKernelGradientMag = kernelGradientsSamples[indexSamples(resolution, r, kh)];
                for(int coord=0 ; coord<3 ; coord++){
                    direction = (pos[coord][particleID]-pos[coord][potNeighborID]) / r;
                    kernelGradients.push_back(direction * currentKernelGradientMag);
                }
            }
        }
    }
}

//*
// Gives the list of the surrounding boxes
void surroundingBoxes(int box, int nBoxesX, int nBoxesY, int nBoxesZ, std::vector<int> &surrBoxes){
    int index_x, index_y, index_z;
    index_x = box/(nBoxesZ*nBoxesY);
    index_y = (box-index_x*nBoxesZ*nBoxesY)/nBoxesZ;
    index_z = box-index_x*nBoxesZ*nBoxesY-index_y*nBoxesZ;

    std::vector<int> tab(6, 1); // Initialized to 1
    std::vector<int> value(9, 0); // Initialized to 0
    value[0]=-1; value[2]=1;
    value[3]=-nBoxesZ; value[5]=nBoxesZ;
    value[6]=-nBoxesZ*nBoxesY; value[8]=nBoxesZ*nBoxesY;

    // Filling of the tab vector.
    if (index_z>0){tab[0]=0;}
    if (index_z<nBoxesZ-1){tab[1]=2;}
    if (index_x>0){tab[4]=0;}
    if (index_x<nBoxesX-1){tab[5]=2;}
    if (index_y>0){tab[2]=0;}
    if (index_y<nBoxesY-1){tab[3]=2;}

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

// Gives the distance to the square between two particles
double distance(std::vector<double> (&pos)[3], int partA, int partB){
    return (pos[0][partA]-pos[0][partB])*(pos[0][partA]-pos[0][partB])
             + (pos[1][partA]-pos[1][partB])*(pos[1][partA]-pos[1][partB])
             + (pos[2][partA]-pos[2][partB])*(pos[2][partA]-pos[2][partB]);

}
