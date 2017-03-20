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
    switch(parameter->viscosityModel)
    {
        case violeauArtificial :
        for (int i = 0; i < neighbors.size(); i++)
        {
            double ux =  currentField->speed[3*particleID]-currentField->speed[3*neighbors[i]];
            double uy =  currentField->speed[3*particleID+1]-currentField->speed[3*neighbors[i]+1];
            double uz =  currentField->speed[3*particleID+2]-currentField->speed[3*neighbors[i]+2];
            double rx =  currentField->pos[3*particleID]-currentField->pos[3*neighbors[i]];
            double ry =  currentField->pos[3*particleID+1]-currentField->pos[3*neighbors[i]+1];
            double rz =  currentField->pos[3*particleID+2]-currentField->pos[3*neighbors[i]+2];
            double Rij_Uij = ux*rx+ry*uy+rz*uz;// moyen plus élégant de faire tout ça ?
            if(Rij_Uij < 0)
            {
                double Rij2 = rx*rx+ry*ry+rz*rz;
                double nu2 = parameter->epsilon*h*h;
                double mu = (h*Rij_Uij)/(Rij2+nu2);
                double rho  = 0.5 * (currentField->density[particleID] + currentField->density[neighbors[i]]);
                viscosity[i] = ( -parameter->alpha*parameter->c*mu + parameter->beta*mu*mu ) / (rho);
            }
            else
            {
                viscosity[i]=0;
            }
        }
        break;

        default :
        viscosity.assign(neighbors.size(),0.0);
        break;
    }

}
