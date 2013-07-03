#ifndef __vtkParallelTree_h__
#define __vtkParallelTree_h__

#include "vtkParallelSMPModule.h" // For export macro
#include "vtkObject.h"

class vtkFunctor;
class vtkLocalData;

//======================================================================================
class VTKPARALLELSMP_EXPORT vtkParallelTree
{
public:
  virtual int TraverseNode( vtkIdType id, int lvl, vtkFunctor* function, vtkLocalData* data ) const = 0;
  virtual void GetTreeSize ( int& max_level, vtkIdType& branching_factor ) const = 0;
};

#endif
// VTK-HeaderTest-Exclude: vtkParallelTree.h
