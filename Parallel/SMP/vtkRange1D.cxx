#include "vtkRange1D.h"
#include "vtkObjectFactory.h"

vtkStandardNewMacro(vtkRange1D);

vtkRange1D::vtkRange1D() : begin(0), end(0)
  {
  tid = 0;
  }

vtkRange1D::~vtkRange1D() {}

void vtkRange1D::PrintSelf(ostream& os, vtkIndent indent)
  {
  this->Superclass::PrintSelf(os,indent);
  }
