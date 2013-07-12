#include "vtkTreeFunctor.h"

//======================================================================================
vtkTreeFunctor::vtkTreeFunctor()
  {
  this->ComputeMasterTID();
  }

//--------------------------------------------------------------------------------
vtkTreeFunctor::~vtkTreeFunctor() { }

//--------------------------------------------------------------------------------
void vtkTreeFunctor::PrintSelf(ostream &os, vtkIndent indent)
  {
  this->Superclass::PrintSelf( os, indent );
  }

//--------------------------------------------------------------------------------
vtkLocalData* vtkTreeFunctor::getLocal(int tid) const
  {
  vtkLocalData* data = vtkLocalData::New();
  return data;
  }
