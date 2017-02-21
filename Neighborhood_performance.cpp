#include "SPH.h"
#include <fstream>
#include <vector>
#include <ctime>

int main(int argc, char *argv[])
{
    //Should be in argv[]
    double o[3] = {0.0,0.0,0.0};
    double L[3] = {20.0,15.0,10.0};
    double s = 1.4;
    double kh = 2.2;
    //End input (end of what should be in argv)

    std::vector<double> pos;
    std::vector<double> values;
    std::vector<int> row;
    std::vector<int> column;

    //Generate cube
    meshcube(o,L,s,pos);

    //Record algorithm performance
    std::clock_t start;
    double duration;

    start = std::clock();
    neighborAllPair(pos, kh, values, row, column);
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout<<"Elapsed time AllPair: "<< duration <<" [s]\n";

    start = std::clock();
    neighborLinkedList(pos,o, L, kh, values, row, column);
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout<<"Elapsed time Linked List: "<< duration <<" [s]\n";



    return 0;
}
