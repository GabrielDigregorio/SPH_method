///**************************************************************************
/// SOURCE: Functions to initialise a field and compute pressure from density.
///**************************************************************************
#include "Main.h"
#include "Physics.h"
/*
*Input:
*- field: field whose speeds will be initialised
*- parameter: pointer the the structure containing parameters
*Decscription:
*Initialise speed from field.
*/
void speedInit(Field* field,Parameter* parameter)
{
	field->speed.assign(3*field->nTotal,0.0); // Initial state is zero speed; other choice could be implemented
  if(field->nMoving != 0){updateMovingSpeed(field,parameter,0.0);}
}

/*
*Input:
*- field: field whose densities will be initialised
*- parameter: pointer the the structure containing parameters
*Decscription:
*Initialise densities from field.
*/
void densityInit(Field* field,Parameter* parameter)
{
	//Parameter withdrawal
	double rho_0 = parameter->densityRef;
	double B = parameter->B;
	double gamma = parameter->gamma;
	double g = parameter->g;
	//std::vector<double> z;
	double zMax = 0.0;

	switch (parameter->densityInitMethod)
	{
		case hydrostatic:

		double H;
		
		for (int j = 0; j < field->nFree; j += 3)
		{
			if (field->pos[j] > zMax)
			{
				zMax = field->pos[j];
			}
			

		}

		for (int i = 0; i < field->nFree; i++)
		{
			H = zMax - field->pos[3 * i + 2];
			double rho = (1 + (1 / B)*rho_0*g*H);
			field->density.push_back(rho_0*pow(rho, 1.0 / gamma));
		}
		for (int k = field->nFree; k < field->nTotal; k++)
		{
			field->density.push_back(parameter->densityRef);
		}
		break;

		case homogeneous:
			field->density.assign(field->nTotal,rho_0);
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
}

/*
*Input:
*- field: field whose pressure will be initialised
*- parameter: pointer the the structure containing parameters
*Decscription:
*Initialise pressure from field.
*/
void pressureInit(Field* field,Parameter* parameter)
{
		field->pressure.resize(field->nTotal);
		pressureComputation(field,parameter);
}

/*
*Input:
*- field: field whose pressure will be updated
*- parameter: pointer the the structure containing parameters
*Decscription:
*Compute pressure from field.
*/
void pressureComputation(Field* field,Parameter* parameter)
{
	//Parameter withdrawal
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

		case perfectGas:
		std::cout << "Perfect gas not yet implemented\n" << std::endl;
		break;

	}
}

/*
*Input:
*- field: field whose masses are initialised
*- parameter: pointer the the structure containing parameters
*Decscription:
*Initialise mass from field.
*/
void massInit(Field* field,Parameter* parameter, std::vector<double> &vol)
{
	for (int i=0; i<field->nTotal; i++)
	{
		double m = field->density[i]*vol[i];
		field->mass.push_back(m);
	}
}
