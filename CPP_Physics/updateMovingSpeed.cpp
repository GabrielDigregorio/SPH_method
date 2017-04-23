#include "Main.h"
#include "Physics.h"
#define M_PI           3.14159265358979323846  /* pi */

double angleComputation(double amplitude,double charactTime,AngleLaw angleLaw, double t)
{
  switch (angleLaw)
  {
    case linearAngle:
    return 2.0*M_PI/charactTime*t;
    break;
    case sineAngle:
    return amplitude*sin(2.0*M_PI/charactTime*t);
    break;
    case exponentialAngle:
    return amplitude*(1.0-exp(-t/charactTime));
    break;
  }
}

double angleDerivativeComputation(double amplitude,double charactTime,AngleLaw angleLaw,double t)
{
  switch (angleLaw)
  {
    case linearAngle:
    return 2.0*M_PI/charactTime;
    break;
    case sineAngle:
    return 2.0*M_PI/charactTime*amplitude*cos(2*M_PI/charactTime*t);
    break;
    case exponentialAngle:
    return amplitude/charactTime*exp(-t/charactTime);
    break;
  }
}
void timeDerivativeQuaternionRotation(double* mD, double angle, double angleDerivative, double* pos)
{
  double s = sin(angle)*angleDerivative;
  double c = cos(angle)*angleDerivative;
  double posTmp[3] = {0,0,0};
  double RDot[3][3] =
  {
    {-s + mD[0]*mD[0]*s, mD[0]*mD[1]*s-mD[2]*c, mD[0]*mD[2]*s+mD[1]*c},
    {mD[0]*mD[1]*s+mD[2]*c, -s + mD[1]*mD[1]*s, mD[1]*mD[2]*s - mD[0]*c},
    {mD[0]*mD[2]*s-mD[1]*c,  mD[1]*mD[2]*s + mD[0]*c, -s + mD[2]*mD[2]*s}
  };
  for(int i=0; i<3;i++)
  {
    for(int j=0;j<3;j++)
    {
      posTmp[i] += RDot[i][j]*pos[j];
    }
  }
  for(int i=0; i<3;i++)
  {
    pos[i] = posTmp[i];
  }
}

void quaternionRotation(double* mD, double angle, double* pos)
{
  double s = sin(angle);
  double c = cos(angle);
  double posTmp[3] = {0,0,0};
  double R[3][3] =
  {
    {c + mD[0]*mD[0]*(1.0-c), mD[0]*mD[1]*(1.0-c)-mD[2]*s, mD[0]*mD[2]*(1.0-c)+mD[1]*s},
    {mD[0]*mD[1]*(1.0-c)+mD[2]*s, c + mD[1]*mD[1]*(1.0-c), mD[1]*mD[2]*(1.0-c) - mD[0]*s},
    {mD[0]*mD[2]*(1.0-c)-mD[1]*s,  mD[1]*mD[2]*(1.0-c) + mD[0]*s, c + mD[2]*mD[2]*(1.0-c)}
  };
  for(int i=0; i<3;i++)
  {
    for(int j=0;j<3;j++)
    {
      posTmp[i] += R[i][j]*pos[j];
    }
  }
  for(int i=0; i<3;i++)
  {
    pos[i] = posTmp[i];
  }
}

/*
 * In: field = structure containing the speed of particules (among others)
 *     parameter = structure containing the parameter usefull to know the movement of the wall
 * Out: Mise à jour des vitesses des parois mobiles
 */
 void updateMovingSpeed(Field* field, Parameter* parameter, double t,double k, int particleID)
 {
  int movingBoundaryID = field->type[particleID]-2;
  double mD[3]= {parameter->movingDirection[0][movingBoundaryID],parameter->movingDirection[1][movingBoundaryID],parameter->movingDirection[2][movingBoundaryID]};
  double rC[3]= {parameter->rotationCenter[0][movingBoundaryID],parameter->rotationCenter[1][movingBoundaryID],parameter->rotationCenter[2][movingBoundaryID]};
  double charactTime=parameter->charactTime[movingBoundaryID];
  double amplitude = parameter->amplitude[movingBoundaryID];
  PosLaw posLaw = (PosLaw) parameter->posLaw[movingBoundaryID];
  AngleLaw angleLaw = (AngleLaw) parameter->angleLaw[movingBoundaryID];

  switch(posLaw)
   {
     case sine:
     {
       // uniform periodic movement in each direction
         // x
         field->speed[0][particleID]=amplitude*mD[0]*2*M_PI/charactTime*cos(2*M_PI/charactTime*(t+k));
         // y
         field->speed[1][particleID]=amplitude*mD[1]*2*M_PI/charactTime*cos(2*M_PI/charactTime*(t+k));
         // z
         field->speed[2][particleID]=amplitude*mD[2]*2*M_PI/charactTime*cos(2*M_PI/charactTime*(t+k));
     }
     break;
     case constant:
     {
       // uniform movement in each direction
         // x
         field->speed[0][particleID]=amplitude*mD[0];
         // y
         field->speed[1][particleID]=amplitude*mD[1];
         // z
         field->speed[2][particleID]=amplitude*mD[2];
     }
     break;
     case exponential:
     // x
     field->speed[0][particleID]=amplitude/charactTime*mD[0]*exp(-(t+k)/charactTime);
     // y
     field->speed[1][particleID]=amplitude/charactTime*mD[1]*exp(-(t+k)/charactTime);
     // z
     field->speed[2][particleID]=amplitude/charactTime*mD[2]*exp(-(t+k)/charactTime);
     break;
     case rotating:
     {
       double pos[3] = {field->pos[0][particleID]-rC[0],field->pos[1][particleID]-rC[1],field->pos[2][particleID]-rC[2]};
       double deltaAngle = angleComputation(amplitude,charactTime,angleLaw,t+k)-angleComputation(amplitude,charactTime,angleLaw,t);
       double angleDerivative = angleDerivativeComputation(amplitude,charactTime,angleLaw,t+k);
       timeDerivativeQuaternionRotation(mD, deltaAngle, angleDerivative, pos);
       for(int i = 0;i<3;i++)
       {
         field->speed[i][particleID] = pos[i];
       }
     }
     break;
   }

 }

 /*
  * In: field = structure containing the speed of particules (among others)
  *     parameter = structure containing the parameter usefull to know the movement of the wall
  * Out: Mise à jour des vitesses des parois mobiles
  */
  void updateMovingPos(Field* field, Parameter* parameter, double t,double k, int particleID)
  {
   int movingBoundaryID = field->type[particleID]-2;
   double mD[3]= {parameter->movingDirection[0][movingBoundaryID],parameter->movingDirection[1][movingBoundaryID],parameter->movingDirection[2][movingBoundaryID]};
   double rC[3]= {parameter->rotationCenter[0][movingBoundaryID],parameter->rotationCenter[1][movingBoundaryID],parameter->rotationCenter[2][movingBoundaryID]};
   double charactTime=parameter->charactTime[movingBoundaryID];
   double amplitude = parameter->amplitude[movingBoundaryID];
   PosLaw posLaw = (PosLaw) parameter->posLaw[movingBoundaryID];
   AngleLaw angleLaw = (AngleLaw) parameter->angleLaw[movingBoundaryID];

   switch(posLaw)
    {
      case sine:
      {
        // uniform periodic movement in each direction
          // x
          field->pos[0][particleID]+=(amplitude*mD[0]*sin(2*M_PI/charactTime*(t+k))-amplitude*mD[0]*sin(2*M_PI/charactTime*(t)));
          // y
          field->pos[1][particleID]+=(amplitude*mD[1]*sin(2*M_PI/charactTime*(t+k))-amplitude*mD[1]*sin(2*M_PI/charactTime*(t)));
          // z
          field->pos[2][particleID]+=(amplitude*mD[2]*sin(2*M_PI/charactTime*(t+k))-amplitude*mD[2]*sin(2*M_PI/charactTime*(t)));
      }
      break;
      case constant:
      {
        // uniform movement in each direction
          // x
          field->pos[0][particleID]+=amplitude*mD[0]*k;
          // y
          field->pos[1][particleID]+=amplitude*mD[1]*k;
          // z
          field->pos[2][particleID]+=amplitude*mD[2]*k;
      }
      break;
      case exponential:
      // x
      field->pos[0][particleID]+=amplitude*mD[0]*(exp(-(t)/charactTime) - exp(-(t+k)/charactTime));
      // y
      field->pos[1][particleID]+=amplitude*mD[1]*(exp(-(t)/charactTime) - exp(-(t+k)/charactTime));
      // z
      field->pos[2][particleID]+=amplitude*mD[2]*(exp(-(t)/charactTime) - exp(-(t+k)/charactTime));
      break;
      case rotating:
      {
        double pos[3] = {field->pos[0][particleID]-rC[0],field->pos[1][particleID]-rC[1],field->pos[2][particleID]-rC[2]};
        double deltaAngle = angleComputation(amplitude,charactTime,angleLaw,t+k)-angleComputation(amplitude,charactTime,angleLaw,t);
        quaternionRotation(mD, deltaAngle, pos);
        for(int i = 0;i<3;i++)
        {
          field->pos[i][particleID] = pos[i] + rC[i];
        }
      }
      break;
    }

  }
