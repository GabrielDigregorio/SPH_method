///**************************************************************************
/// SOURCE: Function to implement the several viscosity equations.
///**************************************************************************
#include "Main.h"
#include "Physics.h"

/*
*Input:
*- particleID: ID of the particle for which equation is computed
*- neighbors: vector of the ID of the particles surrounding the particle for which the equation is computed
*- currentField: field containing information about all particles
*- parameter: user defined parameter stored in a structure
*- viscosity: vector that is filled with the computed viscosities
*Decscription:
*Computes an artificial viscosity and stores it in the vector viscosity
*/
void viscosityComputation(int particleID, std::vector<int>& neighbors, Field* currentField, Parameter* parameter,std::vector<double>& viscosity)
{
    double h = gethFromkh(parameter->kernel ,parameter->kh);
    double maxMu = 0.0;

    switch(parameter->viscosityModel)
    {
        case violeauArtificial :
        for (int i = 0; i < neighbors.size(); i++)
        {
            double ux =  currentField->speed[0][particleID]-currentField->speed[0][neighbors[i]];
            double uy =  currentField->speed[1][particleID]-currentField->speed[1][neighbors[i]];
            double uz =  currentField->speed[2][particleID]-currentField->speed[2][neighbors[i]];
            double rx =  currentField->pos[0][particleID]-currentField->pos[0][neighbors[i]];
            double ry =  currentField->pos[1][particleID]-currentField->pos[1][neighbors[i]];
            double rz =  currentField->pos[2][particleID]-currentField->pos[2][neighbors[i]];
            double Rij_Uij = ux*rx+ry*uy+rz*uz;// moyen plus élégant de faire tout ça ?
            if(Rij_Uij < 0.0)
            {
                double Rij2 = rx*rx+ry*ry+rz*rz;
                double nu2 = parameter->epsilon*h*h;
                double mu = (h*Rij_Uij)/(Rij2+nu2);
                double rho  = 0.5 * (currentField->density[particleID] + currentField->density[neighbors[i]]);

                viscosity[i] = ( -parameter->alpha*parameter->c*mu + parameter->beta*mu*mu ) / (rho);
                if (maxMu < mu){maxMu = mu;}
            }
            else
            {
                viscosity[i]=0.0;
            }
        }
        break;

        default :
        viscosity.assign(neighbors.size(),0.0);
        break;

    }

    // Adaptative Time Step
    switch(parameter->adaptativeTimeStep)
    {
        case yes :
        {
            std::cout << "TIME STEP MUST BE CONSTANT IN THE CURRENT VERSION" << std::endl;

            double t_f  = 0.25 * sqrt(h/parameter->g);
            double t_cv = 0.4  * (h/(parameter->c+0.6*parameter->alpha*parameter->c+0.6*parameter->beta*maxMu));

            if (t_f < t_cv)
                currentField->nextK = t_f ;
            else if(t_cv < t_f)
                currentField->nextK = t_cv ;
        }
        break;
    }
}
