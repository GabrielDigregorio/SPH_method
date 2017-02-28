#include "SPH.hpp"
#include <fstream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <cmath>

int main(int argc, char *argv[])
{
    // no stack geometry
    bool stack = false;

    // open a file to write the time (import into matlab) MUST BE REMOVE LATER
    //ofstream myfile;
    //myfile.open("neighborAnalysis.txt", std::ofstream::out | std::ofstream::app);

    // increment the length of the cube
    //for(int i=1 ; i<=11 ; ++i)
    //{
        //Input parameters
        double s = atof(argv[1]);
        double kh = atof(argv[2]);
        double l = atof(argv[3]);
        double eps = atof(argv[4]);
        std::cout << "\n Parameter list: " << s << ", " << kh << ", " << l << ", " << eps << "\n";

        double o[3] = {0.0,0.0,0.0};
        double L[3] = {l,5,5};

        std::vector<double> pos;
        std::vector<double> valuesNaive;
        std::vector<int> rowNaive;
        std::vector<int> columnNaive;

        std::vector<double> valuesLL;
        std::vector<int> rowLL;
        std::vector<int> columnLL;

        //Generate cube
        meshcube(o,L,s,pos, eps, stack);
        //myfile << (l+1)*(l+1)*(l+1) << " " ;

        double ll[3] = {-L[0]/2, -L[1]/2, -L[2]/2};
        double uu[3] = {L[0]/2, L[1]/2, L[2]/2};

        //Record algorithm performance
        std::clock_t start;
        double duration;

        start = std::clock();
        neighborAllPair(pos, kh, valuesNaive, rowNaive, columnNaive);



        duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
        std::cout<<"Elapsed time AllPair: " << duration <<" [s]\n";
        //myfile << duration << " " ;

        start = std::clock();

        neighborLinkedList(pos, ll, uu, kh, valuesLL, rowLL, columnLL);

        duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
        std::cout<<"Elapsed time Linked List: " << duration <<" [s]\n";
        //myfile << duration << " " << "\n" ;

        std::cout<<"Neighbor pairs for Naive:" << valuesNaive.size() << "\n";
        std::cout<<"Neighbor pairs for Linked-list:" << valuesLL.size() << "\n";

    //}

    std::cout << "\n";
    //myfile.close();

    return 0;
}
