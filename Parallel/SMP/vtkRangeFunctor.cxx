#include "vtkRangeFunctor.h"

//======================================================================================
vtkRangeFunctor::vtkRangeFunctor()
  {
  this->ComputeMasterTID();
  }

//--------------------------------------------------------------------------------
vtkRangeFunctor::~vtkRangeFunctor() { }

//--------------------------------------------------------------------------------
void vtkRangeFunctor::PrintSelf(ostream &os, vtkIndent indent)
  {
  this->Superclass::PrintSelf( os, indent );
  }
