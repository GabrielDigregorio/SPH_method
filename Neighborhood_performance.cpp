#include "SPH.h"
#include <fstream>
#include <vector>
#include <ctime>

int main(int argc, char *argv[])
{
    //Should be in argv[]
    double o[3] = {0.0,0.0,0.0};
    double L[3] = {10.0,10.0,10.0};
    double s = 0.5;
    double kh = 2.5;
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

    std::cout << "All Linked List :\n";
    neighborLinkedList(pos,o, L, kh, values, row, column);
    std::cout << "All pair :\n";
    //neighborAllPair(pos, kh, values, row, column);

    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

    std::cout<<"Elapsed time: "<< duration <<" [s]\n";

    return 0;
}
