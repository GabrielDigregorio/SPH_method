#include "SPH.h"
#include <fstream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <cmath>
using namespace std;

int main(int argc, char *argv[])
{
    // open a file to write the time (import into matlab) MUST BE REMOVE LATER
    ofstream myfile;
    myfile.open("neighborAnalysis.txt", std::ofstream::out | std::ofstream::app);

for(int i=1 ; i<=11 ; ++i)
{
    //Input parameters
    double s = atof(argv[1]);
    double kh = atof(argv[2]);
    double l = atof(argv[3]) + i ;
    std::cout << "\n Parameter list: " << s << ", " << kh << ", " << l << "\n";

    double o[3] = {0.0,0.0,0.0};
    double L[3] = {l,l,l};

    std::vector<double> pos;
    std::vector<double> values;
    std::vector<int> row;
    std::vector<int> column;

    //Generate cube
    meshcube(o,L,s,pos, 0);
    myfile << (l+1)*(l+1)*(l+1) << " " ;

    //Record algorithm performance
    std::clock_t start;
    double duration;

    start = std::clock();
    neighborAllPair(pos, kh, values, row, column);
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout<<"Elapsed time AllPair: "<< duration <<" [s]\n";
    myfile << duration << " " ;

    start = std::clock();
    double ll[3] = {-L[0]/2, -L[1]/2, -L[2]/2};
    double uu[3] = {L[0]/2, L[1]/2, L[2]/2};
    neighborLinkedList(pos, ll, uu, kh, values, row, column);
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout<<"Elapsed time Linked List: "<< duration <<" [s]\n";
    myfile << duration << " " << "\n" ;
}


    std::cout << "\n";
    myfile.close();
    return 0;
}
