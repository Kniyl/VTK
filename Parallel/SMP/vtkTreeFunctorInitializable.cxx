#include "vtkTreeFunctorInitializable.h"

extern int vtkSMPInternalGetNumberOfThreads();
extern int vtkSMPInternalGetTid();

//======================================================================================
vtkTreeFunctorInitializable::vtkTreeFunctorInitializable() :
    vtkTreeFunctor(),
    IsInitialized(vtkSMPInternalGetNumberOfThreads(), 0)
  {
  }

//--------------------------------------------------------------------------------
vtkTreeFunctorInitializable::~vtkTreeFunctorInitializable()
  {
  IsInitialized.clear();
  }

//--------------------------------------------------------------------------------
bool vtkTreeFunctorInitializable::ShouldInitialize(int tid) const
  {
  return !IsInitialized[tid];
  }

//--------------------------------------------------------------------------------
void vtkTreeFunctorInitializable::Initialized(int tid) const
  {
  IsInitialized[tid] = 1;
  }

//--------------------------------------------------------------------------------
void vtkTreeFunctorInitializable::PrintSelf(ostream &os, vtkIndent indent)
  {
  this->Superclass::PrintSelf( os, indent );
  os << indent << "Is initialized: " << endl;
  for ( vtkstd::vector<vtkIdType>::size_type i = 0; i < IsInitialized.size(); ++i )
    {
    os << indent.GetNextIndent() << "Id " << i << ": ";
    if ( IsInitialized[i] )
      os << "true";
    else
      os << "false";
    os << endl;
    }
  }
