#include "vtkFunctor.h"

//======================================================================================
vtkFunctor::vtkFunctor()
  {
  this->ComputeMasterTID();
  }

//--------------------------------------------------------------------------------
vtkFunctor::~vtkFunctor() { }

//--------------------------------------------------------------------------------
void vtkFunctor::PrintSelf(ostream &os, vtkIndent indent)
  {
  this->Superclass::PrintSelf( os, indent );
  }

//--------------------------------------------------------------------------------
vtkLocalData* vtkFunctor::getLocal(int tid) const
  {
  vtkLocalData* data = vtkLocalData::New();
  return data;
  }
