#include "Main.h"
#include "Physics.h"

void speedInit(Field* field,Parameter* parameter)
{
	std::cout << "----BEGIN speed initialitation---- \n \n";
	// On démarre du repos, d'autres choix pourraient être fait.
	field->speed.assign(field->nTotal,0.0);

  // Initialisation des moving boundaries
  if(field->nMoving != 0){updateMovingSpeed(field,parameter,0.0);}

	std::cout << "----END speed initialitation---- \n \n";
}
//-----------------------------------------------------------------------------------------
/*
* In: field = structure containing the position of particules (among others)
*     parameter = structure containing the parameter usefull to compute the densities
* Out: Initialisation of the density for the field
*/
void densityInit(Field* field,Parameter* parameter)
{
	std::cout << "----BEGIN density initialitation---- \n \n";
	//Récupération des paramètres
	double rho_0 = parameter->densityRef;
	double B = parameter->B;
	double gamma = parameter->gamma;
	double g = parameter->g;
	switch (parameter->densityInitMethod)
	{
		case hydrostatic:

		double H;
		// Find height of free surface.
		double zMax;
		zMax = *(std::max_element(&field->pos[0],&field->pos[field->nFree+1]));

		for (int i = 0; i < field->nTotal; i++)
		{
			H = zMax - field->pos[3*i+2];
			double rho = rho_0*(1 + (1 / B)*rho_0*g*H);

			field->density[i] = pow(rho, 1.0 / gamma);
		}
		break;
		case homogeneous:
		for (int i = 0; i < field->nTotal; i++)
		{
			field->density[i] = rho_0;
		}

		break;
		// A voir si on considère le cas d'un gas pft, nécéssite des parametres supplémentaires(T°,M)
		/*
		if (parameter->stateEquationMethod == perfectGas)
		{
		for (int i=2; i<field->pos.size(); i=i+3)
		{
		H = z_max - field->pos[i];
		double rho = rho_0*(1 + (Mach/R*Temperature)*rho_0*g*H);
		field->density.push_back(rho)
	}

}
*/
}
std::cout << "----END density initialitation---- \n \n";
}

//-----------------------------------------------------------------------------------------
/*
* In: field = structure containing the density of particules (among others)
*     parameter = structure containing the parameter usefull to compute the pressures
* Out: Computation of the pressure for the field
*/
void pressureComputation(Field* field,Parameter* parameter)
{

	std::cout << "----BEGIN pressure computation---- \n \n";
	//Récupération des paramètres
	double rho_0 = parameter->densityRef;
	double B = parameter->B;
	double gamma = parameter->gamma;
	double g = parameter->g;

	switch (parameter->stateEquationMethod)
	{
		case quasiIncompressible:


		for (int i=0; i<field->nTotal; i++)
		{
			double rho = field->density[i];
			double p = B*(pow(rho/rho_0,gamma)-1);

			field->pressure[i] = p;
		}
		break;
	}

	std::cout << "----END pressure computation---- \n \n";
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

	std::cout << "----BEGIN mass initialitation---- \n \n";
	for (int i=0; i<field->nTotal; i++)
	{
		double V = field->s[i]*field->s[i]*field->s[i];
		double m = field->density[i]*V;

		field->mass[i] = m;
	}

	std::cout << "----END mass initialitation---- \n \n";

}
