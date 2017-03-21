#include "Main.h"
#include "Interface.h"

// Clear the boxes content
void boxClear(std::vector<std::vector<int> > &boxes)
{
    for(int i=0 ; i<boxes.size() ; i++)
        boxes[i].clear();
}

/*
*Input:
*- sourceField: pointer to the field to be copied
*- copiedField: pointer to an empty field
*Decscription:
*Copy unvariant datas of a field structure and give the right size to time varying datas
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

  copiedField->pos.resize(3*nTotal);
  copiedField->speed.resize(3*nTotal);
  copiedField->pressure.resize(nTotal);
  copiedField->density.resize(nTotal);

  // Copying fixed positions
  for(int i = 3*sourceField->nFree; i < 3*(sourceField->nFree + sourceField->nFixed);i++)
  {
    copiedField->pos[i] = sourceField->pos[i];
  }
}

/*
*Input:
*- hopField/cornField: fields to swap
*Description:
*Exchange the content of the two fields.
*/
void swapField(Field** hopField, Field** cornField)
{
  Field *tmpField;
  tmpField = *hopField;
  *hopField = *cornField;
  *cornField = tmpField;
}
