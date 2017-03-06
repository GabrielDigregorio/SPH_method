#include "../Headers/SPH.hpp"

//-----------------------------------------------------------------------------------------
/*
 * In: field = structure containing the position of particules (among others)
 *     parameter = structure containing the parameter usefull to compute the densities
 * Out: Initialisation of the density for the field
 */
void densityInit(Field* field,Parameter* parameter)
{
    //Récupération des paramètres
    double rho_0 = parameter->densityRef;
    double B = parameter->B;
    double gamma = parameter->gamma;
    double g = parameter->g; 

    double H;
    double zMax = field->u[2];
   
    //Initialistation of density for each particles(Case of quasiIncompressible)
    for (int i=2; i<field->posFree.size(); i=i+3)
    {
        H = zMax - field->posFree[i];
        rho = rho_0*(1 + (1/B)*rho_0*g*H);

        field->densityFree.push_back(pow(rho,1/gamma));
    }

    for (int j=2; j<field->posMoving.size(); j=j+3)
    {
        H = zMax - field->posMoving[j];
        rho = rho_0*(1 + (1/B)*rho_0*g*H);

        field->densityMoving.push_back(pow(rho,1/gamma));
    }

    for (int k=2; k<field->posFixed.size(); k=k+3)
    {
        H = zMax - field->posFixed[k];
        rho = rho_0*(1 + (1/B)*rho_0*g*H);

        field->densityFixed.push_back(pow(rho,1/gamma));
    }


    

 // A voir si on considère le cas d'un gas pft, nécéssite des parametres supplémentaires(T°,M)   
/*
    if (parameter->stateEquationMethod == "perfectGas")
    {
        for (int i=0; i<field->posFree.size()/3; i++)
        {
            double rho = rho_0*(1 + (Mach/R*Temperature)*rho_0*g*H_free[i]); 
            field->densityFree.push_back(rho)
        } 

        for (int j=0; j<field->posMoving.size()/3; j++)
        {
            double rho = rho_0*(1 + (Mach/R*Temperature)*rho_0*g*H_moving[i]); 
            field->densityFree.push_back(rho)
        } 

        for (int k=0; k<field->posFixed.size()/3; k++)
        {
            double rho = rho_0*(1 + (Mach/R*Temperature)*rho_0*g*H_fixed[i]); 
            field->densityFree.push_back(rho)
        } 
    }
*/ 

}

//-----------------------------------------------------------------------------------------
/*
 * In: field = structure containing the density of particules (among others)
 *     parameter = structure containing the parameter usefull to compute the pressures
 * Out: Computation of the pressure for the field
 */
void pressureComputation(Field* field,Parameter* parameter)
{
    //Récupération des paramètres
    double rho_0 = parameter->densityRef;
    double B = parameter->B;
    double gamma = parameter->gamma;
    double g = parameter->g;

    //Cas d'un liquide quasi_incompressible 

    for (int i=0; i<field->posFree/3; i++)
    {
        double rho = field->densityFree[i];
        double p = B*(pow(rho/rho_0,gamma)-1);

        field->pressureFree.push_back(p);
    }

    for (int j=0; j<field->posMoving/3; j++)
    {
        double rho = field->densityMoving[j];
        double p = B*(pow(rho/rho_0,gamma)-1);

        field->pressureMoving.push_back(p);
    }

    for (int k=0; k<field->posFixed/3; k++)
    {
        double rho = field->densityFixed[k];
        double p = B*(pow(rho/rho_0,gamma)-1);

        field->pressureFixed.push_back(p);
    }

    // + cas gaz parfait ?
}


//-----------------------------------------------------------------------------------------
/*
 * In: field = structure containing the density of particules (among others)
 *     parameter = structure containing the parameter s with is needed to compute the masses
 * Out: Computation of the masses for the field
 */
void massInit(Field* field,Parameter* parameter)
{
    for (int i=0; i<field->posFree/3; i++)
    {
        double V = pow(field->sFree[0],3);
        double m = field->densityFree[i]*V;
        
        field->massFree.push_back(m);
    }

    for (int j=0; j<field->posMoving/3; j++)
    {
        double V = pow(field->sMoving[0],3);
        double m = field->densityMoving[j]*V;
        
        field->massMoving.push_back(m);
    }

    for (int k=0; k<field->posFixed/3; k++)
    {
        double V = pow(field->sFixed[0],3);
        double m = field->densityFixed[k]*V;
        
        field->massFixed.push_back(m);
    }

}
