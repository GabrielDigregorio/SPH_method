#ifndef PLAYGROUND_H
#define PLAYGROUND_H

// Structure qui doit être remplie lors de la lecture du fichier de paramètre
/*
 * kh = smothing length
 * k = time step
 * T = simulation time
 * densityRef = density of the fluid at atmospheric pressure
 * l & u = lower and upper limit of the domain
 * B & gamma = fluid constants
 * g = gravity
 * writeInteval = time interval between two outputs file are generated
 * integrationMethod = euler ou RK2
 * densityInitMethod = hydrosatic, etc.
 * stateEquationMethod = quasiIncompressible, perfectGas, etc.
 * massInitMethod = violeau2012 (all particles have same volumes), etc.
 * speedLaw = To be determined, will dictate the behaviour of moving boundaries
*/

struct Parameter {
    double kh, k, T, densityRef, B, gamma, g, writeInterval;
    double l[3];
    double u[3];
    std::string integrationMethod, densityInitMethod, stateEquationMethod, massInitMethod, speedLaw;
};

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

        // Get parameters of the fluid
        Parameter* GetParam();

        // Destructor: delete all vectors in DATA
        ~Playground();

    private:

        // Print on terminal
        bool screen;

        // Parameters of the fluid
        std::vector<double> param;

        // Method used for the solver (Euler or RungeKutta, ...)
        std::vector< std::string > method;

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