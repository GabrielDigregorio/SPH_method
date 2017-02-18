#include "sph.h"



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
void neighborLinkedList (double p[3], double h, std::vector<double> &pos)
{

}

// tree search algorithm
void neighborTree (double p[3], double h, std::vector<double> &pos)
{

}
