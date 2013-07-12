#include "vtkRangeFunctorInitializable.h"

extern int vtkSMPInternalGetNumberOfThreads();
extern int vtkSMPInternalGetTid();

//======================================================================================
vtkRangeFunctorInitializable::vtkRangeFunctorInitializable() :
    vtkRangeFunctor(),
    IsInitialized(vtkSMPInternalGetNumberOfThreads(), 0)
  {
  }

//--------------------------------------------------------------------------------
vtkRangeFunctorInitializable::~vtkRangeFunctorInitializable()
  {
  IsInitialized.clear();
  }

//--------------------------------------------------------------------------------
bool vtkRangeFunctorInitializable::ShouldInitialize(int tid) const
  {
  return !IsInitialized[tid];
  }

//--------------------------------------------------------------------------------
void vtkRangeFunctorInitializable::Initialized(int tid) const
  {
  IsInitialized[tid] = 1;
  }

//--------------------------------------------------------------------------------
void vtkRangeFunctorInitializable::PrintSelf(ostream &os, vtkIndent indent)
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
