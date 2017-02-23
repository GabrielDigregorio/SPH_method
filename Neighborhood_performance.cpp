#include "SPH.h"
#include <fstream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <cmath>

int main(int argc, char *argv[])
{
    //Input parameters
    double s = atof(argv[1]);
    double kh = atof(argv[2]);
    double l = atof(argv[3]);
    std::cout << "\n Parameter list: " << s << ", " << kh << ", " << l << "\n";

    double o[3] = {0.0,0.0,0.0};
    double L[3] = {l,l,l};

    std::vector<double> pos;
    std::vector<double> values;
    std::vector<int> row;
    std::vector<int> column;

    //Generate cube
    meshcube(o,L,s,pos, 0);

    //Record algorithm performance
    std::clock_t start;
    double duration;

    start = std::clock();
    neighborAllPair(pos, kh, values, row, column);
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout<<"Elapsed time AllPair: "<< duration <<" [s]\n";

    start = std::clock();
    double l[3] = {-L[0]/2, -L[1]/2, -L[2]/2};
    double u[3] = {L[0]/2, L[1]/2, L[2]/2};
    neighborLinkedList(pos, l, u, kh, values, row, column);
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout<<"Elapsed time Linked List: "<< duration <<" [s]\n";


    std::cout << "\n";
    return 0;
}
