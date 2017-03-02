#ifndef PLAYGROUND_H
#define PLAYGROUND_H

// Class Playground
class Playground
{
    public:

    Playground(); // constructor
    void fillVector(double *A, double *B, double *C, int i);
    void ReadPlayground(const char *filename);
    void GeneratePlayground( std::vector<double> &posFree, std::vector<double> &posMoving, std::vector<double> &posFixed);

    private:

    // 2D Vector: [cube cylinder sphere ; part1 part2 part3 part4 ...]. Store geometry
    std::vector< std::vector<int> > geometry;

    // 3D Vector: [Free Moving Fixed ; Param Coord Dimen ; part1 part2 part3 part4 ...]. Store all parameters
    std::vector< std::vector< std::vector<double> > > DATA;
};

#endif 