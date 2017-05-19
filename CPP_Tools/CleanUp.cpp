#include "Main.h"
#include "Interface.h"

// Clear the boxes content
void boxClear(std::vector<std::vector<int> > &boxes)
{
  for(int i=0 ; i<boxes.size() ; i++)
  boxes[i].clear();
}

// Get the time (or not)
std::clock_t getTime(){
  return 0.0;
}



/*
*Input:
*- sourceField: pointer to the field to be copied
*- copiedField: pointer to an empty field
*Decscription:
*Copy ALL (!!) data of a field structure and give the right size to time varying data
*/
void copyField(Field *sourceField,Field *copiedField)
{

  int nTotal = sourceField->nTotal;
  for (int i = 0; i < 3; i++)
  {
    copiedField->l[i] = sourceField->l[i];
    copiedField->u[i] = sourceField->u[i];
  }
  copiedField->nFree = sourceField->nFree;
  copiedField->nFixed = sourceField->nFixed;
  copiedField->nMoving = sourceField->nMoving;
  copiedField->nTotal = nTotal;

  copiedField->mass = sourceField->mass;
  copiedField->type = sourceField->type;
  copiedField->pressure = sourceField->pressure;
  copiedField->density = sourceField->density;

  for(int j=0 ; j<3 ; j++){
    copiedField->pos[j] = sourceField->pos[j];//.resize(nTotal);//
    copiedField->speed[j] = sourceField->speed[j];//.resize(nTotal);
  }
}

/*
*Input:
*- hopField/cornField: fields to swap
*Description:
*Exchange the content of the two fields because rotation of the cultures increases the fertility of fields.
*/
void swapField(Field** hopField, Field** cornField)
{
  Field *tmpField;
  tmpField = *hopField;
  *hopField = *cornField;
  *cornField = tmpField;
}
