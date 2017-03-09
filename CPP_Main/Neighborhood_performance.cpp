#include "Main.h"
#include "Interface.h"
#include "Physics.h"
#include "Tools.h"
#include <ctime>

int main(int argc, char *argv[])
{
    // no stack geometry
    bool stack = false;

    //Input parameters
    if(argc != 5){
        std::cout << "Invalid input parameters. Must be: s kh l eps\n";
        return EXIT_FAILURE;
    }
    double s = atof(argv[1]);
    double kh = atof(argv[2]);
    double l = atof(argv[3]);
    double eps = atof(argv[4]);
    std::cout << "\n Parameter list: " << s << ", " << kh << ", " << l << ", " << eps << "\n";

    double o[3] = {0.0,0.0,0.0};
    double L[3] = {l,4,5};

    std::vector<double> pos;
    std::vector<double> valuesNaive;
    std::vector<int> rowNaive;
    std::vector<int> columnNaive;

    std::vector<double> valuesLL;
    std::vector<int> rowLL;
    std::vector<int> columnLL;

    //Generate cube
    meshcube(o, L, s, pos, eps);

    double ll[3] = {-L[0]/2, -L[1]/2, -L[2]/2};
    double uu[3] = {L[0]/2, L[1]/2, L[2]/2};


    //Record algorithm performance
    // ALL PAIRS
    std::clock_t start;
    double duration;

    start = std::clock();

    neighborAllPair(pos, kh, valuesNaive, rowNaive, columnNaive);

    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout<<"Elapsed time AllPair: " << duration <<" [s]\n";

    // LINKED LIST
    start = std::clock();
    //for(int i=0; i<10 ; i++)
    neighborLinkedList(pos, ll, uu, kh, valuesLL, rowLL, columnLL);

    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout<<"Elapsed time Linked List: " << duration <<" [s]\n";


    // SPLITTED NEIGHBORS
    start = std::clock();
    //timeIntegration(pos, ll, uu, kh);
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout<<"Elapsed time Splitted: " << duration <<" [s]\n";


    // Comparison of the neighbors
    /*
    ...
    */




    std::cout<<"Neighbor pairs for Naive:" << valuesNaive.size() << "\n";
    std::cout<<"Neighbor pairs for Linked-list:" << valuesLL.size() << "\n";
    //std::cout<<"Neighbor pairs for Splitted:" << nPairs << "\n"; OLD FEATURE

    std::cout << "\n";

    return 0;
}
