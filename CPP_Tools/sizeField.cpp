#include "Main.h"
#include "Interface.h"

void sizeField(Field *field, int nTotal)
{
  field->pos.resize(3*nTotal);
  field->speed.resize(3*nTotal);
  field->mass.resize(nTotal);
  field->pressure.resize(nTotal);
  field->density.resize(nTotal);
  field->s.resize(nTotal);
}
