#include "vtkLocalData.h"
#include "vtkObjectFactory.h"

vtkStandardNewMacro(vtkLocalData);

vtkLocalData::vtkLocalData() {}
vtkLocalData::~vtkLocalData() {}

void vtkLocalData::PrintSelf(ostream& os, vtkIndent indent)
  {
  this->Superclass::PrintSelf(os,indent);
  }
