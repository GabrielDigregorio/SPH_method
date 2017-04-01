#include "Main.h"
#include "Physics.h"
#define M_PI           3.14159265358979323846  /* pi */
/*
 * In: field = structure containing the speed of particules (among others)
 *     parameter = structure containing the parameter usefull to know the movement of the wall
 * Out: Mise à jour des vitesses des parois mobiles
 */

void updateMovingSpeed(Field* field, Parameter* parameter, double t)
{
    int start;
    //int end;
    start = field->nFixed + field->nFree;
    //end = field->nTotal;

    int number_block_diff=field->info_block.size()/6;
    int block_number=0;
    int counter_1=0;
    int counter_2=0;
    while(block_number<number_block_diff)
    {
      int c=field->info_block[counter_2];
      if(c==1)
      {
        counter_2=counter_2+1;
        int nx=field->info_block[counter_2];
        counter_2=counter_2+1;
        int ny=field->info_block[counter_2];
        counter_2=counter_2+1;
        int nz=field->info_block[counter_2];
        counter_2=counter_2+1;
        double angle=field->info_block[counter_2];
        counter_2=counter_2+1;
        int direction=field->info_block[counter_2];
        counter_2=counter_2+1;

        int speedLaw=field->info_moving[counter_1];
        counter_1=counter_1+1;
        double movx=field->info_moving[counter_1];
        counter_1=counter_1+1;
        double movy=field->info_moving[counter_1];
        counter_1=counter_1+1;
        double movz=field->info_moving[counter_1];
        counter_1=counter_1+1;

        int p=start;
        int lm=1;
        int cnt_p=1;

        while(lm<=block_number)
        {
          int nx_p=field->info_block[cnt_p];
          int ny_p=field->info_block[cnt_p+1];
          int nz_p=field->info_block[cnt_p+2];
          cnt_p=cnt_p+6;
          p=p+nx_p*nz_p*ny_p;
          lm=lm+1;
        }
        block_number=block_number+1;

    switch(speedLaw)
    {
      case sine:
      {
      // uniform periodic movement in each direction
      int end=p+nx*ny*nz;
      for(int i=p ; i<end ; i++)
      {
        // x
        field->speed[0][i]=movx*sin(parameter->charactTime*t);
        // y
        field->speed[1][i]=movy*sin(parameter->charactTime*t);
        // z
        field->speed[2][i]=movz*sin(parameter->charactTime*t);
      }
      }
      break;
      case constant:
      {
      int end=p+nx*ny*nz;
      for(int i=p ; i<end ; i++)
      {
        // x
        field->speed[0][i]=movx;
        // y
        field->speed[1][i]=movy;
        // z
        field->speed[2][i]=movz;
      }
      }
      break;
      case exponential:
      std::cout << "Exponential moving not yet implemented\n"; // ATTENTION !!!
      break;
      case level_arm:
      {

        double frequency=2;

            for(int k=0; k<nz; ++k)
              {

                for(int j=0; j<ny; ++j)
                {

                    for(int i=0; i<nx; ++i)
                    {

                      double teta=angle/180.0*M_PI*cos(parameter->charactTime*t);
                      double teta_dot=-angle/180.0*M_PI*frequency*sin(parameter->charactTime*t);
                      // ne fonctionne que dans l'approximation des petits angles, sinon il faut utiliser la forme de jacobi
                      // pour le cas du pendule
                      if(direction==2)
                      {
                      // x
                      field->speed[0][i]=k/10.0*cos(teta)*teta_dot;
                      // y
                      field->speed[1][i]=0.0;
                      // z
                      field->speed[2][i]=-k/10.0*sin(teta)*teta_dot;
                      }
                      else if(direction==1)
                      {
                      // x
                      field->speed[0][i]=0.0;
                      // y
                      field->speed[1][i]=k/10.0*cos(teta)*teta_dot;
                      // z
                      field->speed[2][i]=-k/10.0*sin(teta)*teta_dot;
                      }
                      else if(direction==3)
                      {
                      // x
                      field->speed[0][i]=-j/10.0*cos(teta)*teta_dot;
                      // y
                      field->speed[1][i]=j/10.0*sin(teta)*teta_dot;
                      // z
                      field->speed[2][i]=0.0;
                      }
                      else //if (direction==0)
                      {
                        // x
                      field->speed[0][i]=0.0;
                      // y
                      field->speed[1][i]=0.0;
                      // z
                      field->speed[2][i]=0.0;
                      }

                      p=p+1;
                    }
                }
              }


      }// end case
      break;
      }
     }
     else
     {
       counter_1=counter_1+4;
       counter_2=counter_2+6;
       block_number=block_number+1;
     }
    }
  }
