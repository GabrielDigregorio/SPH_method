///**************************************************************************
/// SOURCE: Check the coherence of input values.
///**************************************************************************
#include "Main.h"
#include "Interface.h"

Error consistency(Parameter* parameter, Field* field)
{
  int cntError = 0;
  if( (field->u[0]<field->l[0]) || (field->u[1]<field->l[1]) || (field->u[2]<field->l[2]) )
  {
    std::cout << "Lower and higher domain dimension are not consistent.\n" << std::endl;
    cntError++;
  }

  int cntOutOfDomain = 0;
  for(int i = 0; i < field->nTotal; i++)
  for (int j=0;j<3;j++)
  {
    {
      if(field->pos[j][i] > field->u[j])
      {
        cntOutOfDomain++;
        break;
      }
      else if (field->pos[j][i] < field->l[j])
      {
        cntOutOfDomain++;
        break;
      }
    }
  }
  if(cntOutOfDomain != 0)
  {
    std::cout << cntOutOfDomain << " particles out of the domain.\n" << std::endl;
    cntError++;
  }

  if(parameter->kh <= 0.0)
  {
    std::cout << "Invalid kh.\n" << std::endl;
    cntError++;
  }

  if( (parameter->k <= 0.0) || (parameter->k > parameter->T) )
  {
    std::cout << "Invalid timestep.\n" << std::endl;
    cntError++;
  }

  if(parameter->densityRef <= 0.0)
  {
    std::cout << "Invalid referecne density.\n" << std::endl;
    cntError++;
  }


  if(parameter->B < 0.0)
  {
    std::cout << "Invalid B.\n" << std::endl;
    cntError++;
  }
  if(parameter->gamma < 0.0)
  {
    std::cout << "Invalid gamma.\n" << std::endl;
    cntError++;
  }
  if(parameter->writeInterval < parameter->k)
  {
    std::cout << "Invalid writeInterval.\n" << std::endl;
    cntError++;
  }
  if(parameter->c < 0.0)
  {
    std::cout << "Invalid speed of sound.\n" << std::endl;
    cntError++;
  }
  if(parameter->alpha < 0.0)
  {
    std::cout << "Invalid alpha.\n" << std::endl;
    cntError++;
  }
  if(parameter->beta < 0.0)
  {
    std::cout << "Invalid beta.\n" << std::endl;
    cntError++;
  }
  if(parameter->epsilon < 0.0)
  {
    std::cout << "Invalid epsilon.\n" << std::endl;
    cntError++;
  }
  if(parameter->temperature < 0.0)
  {
    std::cout << "Invalid temperature.\n" << std::endl;
    cntError++;
  }
  if(parameter->molarMass < 0.0)
  {
    std::cout << "Invalid molarMass.\n" << std::endl;
    cntError++;
  }
  if( (parameter->theta <= 0.0) || (parameter->theta > 1.0) )
  {
    std::cout << "Invalid theta.\n" << std::endl;
    cntError++;
  }
  if (cntError != 0)
  {
    return consistencyError;
  }
  else
  {
    return noError;
  }
}
