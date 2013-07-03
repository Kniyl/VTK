#include "vtkFunctorInitializable.h"

extern int vtkSMPInternalGetNumberOfThreads();
extern int vtkSMPInternalGetTid();

//======================================================================================
vtkFunctorInitializable::vtkFunctorInitializable() :
    vtkFunctor(),
    IsInitialized(vtkSMPInternalGetNumberOfThreads(), 0)
  {
  }

//--------------------------------------------------------------------------------
vtkFunctorInitializable::~vtkFunctorInitializable()
  {
  IsInitialized.clear();
  }

//--------------------------------------------------------------------------------
bool vtkFunctorInitializable::ShouldInitialize(int tid) const
  {
  return !IsInitialized[tid];
  }

//--------------------------------------------------------------------------------
void vtkFunctorInitializable::Initialized(int tid) const
  {
  IsInitialized[tid] = 1;
  }

//--------------------------------------------------------------------------------
void vtkFunctorInitializable::PrintSelf(ostream &os, vtkIndent indent)
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
