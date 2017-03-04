#ifndef PLAYGROUND_H
#define PLAYGROUND_H

// Class Playground
class Playground
{
    public:

        // Constructor: initialize 2D vector and 3D vector, and "screen" varliable
        Playground(bool s);

        // Read the entire file .kzr and store in DATA
        void ReadPlayground(const char *filename);

        // Generate paricles in each geometries accordingly to DATA
        void GeneratePlayground( std::vector<double> &posFree, std::vector<double> &posMoving, std::vector<double> &posFixed);

        // Get the Domain of the playground: dom=false for lower coordinate, dom=true for upper coordinate
        std::vector<double> GetDomain(bool dom);

        // Get parameters of the fluid
        std::vector<double> GetParam();

        // Get method for the solver (Euler or RungeKutta)
        std::string GetMethod();

        // Destructor: delete all vectors in DATA
        ~Playground();

    private:

        // Print on terminal
        bool screen;

        // Parameters of the fluid
        std::vector<double> param;

        // Method used for the solver (Euler or RungeKutta)
        std::string method;

        // Domain of the Playground
        std::vector<double> l{0,0,0}, u{0,0,0};

        // 2D Vector: [cube cylinder sphere ; part1 part2 part3 part4 ...]. Store geometry
        std::vector< std::vector<int> > geometry;

        // 3D Vector: [Free Moving Fixed ; Param Coord Dimen ; part1 part2 part3 part4 ...]. Store all parameters
        std::vector< std::vector< std::vector<double> > > DATA;
        
        // Number of data to read per geometry (3 parameters, 3 coordinates, 3 dimensions)
        int nbr_data = 9;

        // Currently read data
        std::vector<double> data{0,0,0,0,0,0,0,0,0};

        // Fill the 3D matrix DATA with data (current reading)
        void fillVector(std::vector<double> data, int i);
};

#endif 