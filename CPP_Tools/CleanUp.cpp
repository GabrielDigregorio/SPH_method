#include "Main.h"
#include "Interface.h"

// Clear Field structure
void cleanField(Field* field)
{
    field->pos.clear();
    field->speed.clear();
    field->density.clear();
    field->pressure.clear();
    field->mass.clear();
    field->s.clear();
    delete field;
}

// Clear Parameter structure
void cleanParameter(Parameter* parameter)
{
    delete parameter;
}

// Clear the boxes content
void boxClear(std::vector<std::vector<int> > &boxes){
    for(int i=0 ; i<boxes.size() ; i++)
        boxes[i].clear();
}

// Allocate memory
void sizeField(Field *field, int nTotal)
{
  field->pos.resize(3*nTotal);
  field->speed.resize(3*nTotal);
  field->mass.resize(nTotal);
  field->pressure.resize(nTotal);
  field->density.resize(nTotal);
  field->s.resize(nTotal);
}