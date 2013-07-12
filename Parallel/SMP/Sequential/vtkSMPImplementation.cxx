#include "vtkParallelOperators.h"
#include "vtkTreeFunctor.h"
#include "vtkTreeFunctorInitializable.h"
#include "vtkRangeFunctor.h"
#include "vtkRangeFunctorInitializable.h"
#include "vtkRange1D.h"
#include "vtkParallelTree.h"
#include "vtkTask.h"
#include "vtkMergeDataSets.h"

void sequential_traverse( vtkIdType index, int lvl, vtkIdType BranchingFactor, const vtkParallelTree* Tree, vtkTreeFunctor* func, vtkLocalData* data )
  {
  if ( Tree->TraverseNode( index, lvl, func, data ) )
    {
    for ( vtkIdType i = index * BranchingFactor + 1, j = 0; j < BranchingFactor; ++i, ++j )
      {
      sequential_traverse( i, lvl + 1, BranchingFactor, Tree, func, data );
      }
    }
  }

int vtkSMPInternalGetTid()
{
  return 0;
}

int vtkSMPInternalGetNumberOfThreads()
{
  return 1;
}

void vtkParallelOperators::ForEach( vtkIdType first, vtkIdType last, const vtkRangeFunctor* op, int grain )
{
  vtkRange1D* range = vtkRange1D::New();
  range->Setup(first,last,0);
  (*op)( range );
  range->Delete();
}

void vtkParallelOperators::ForEach( vtkIdType first, vtkIdType last, const vtkRangeFunctorInitializable* f, int grain )
{
  if ( f->ShouldInitialize(0) )
    f->Init(0);
  vtkRange1D* range = vtkRange1D::New();
  range->Setup(first,last,0);
  (*f)( range );
  range->Delete();
}

void vtkParallelOperators::Traverse( const vtkParallelTree* Tree, vtkTreeFunctor* func )
{
  int lvl;
  vtkIdType bf;
  Tree->GetTreeSize(lvl,bf);
  vtkLocalData* data = func->getLocal(0);
  sequential_traverse( 0, 0, bf, Tree, func, data );
  data->Delete();
}

void vtkMergeDataSets::Parallel(
    const vtkTask* function,
    vtkThreadLocal<vtkSMPMergePoints>::iterator data)
{
  if ( !this->MasterThreadPopulatedOutput )
    function->Execute( data[0] );
}

void vtkMergeDataSets::Parallel(
   const vtkTask* function,
   vtkThreadLocal<vtkIdList>::iterator data1,
   vtkThreadLocal<vtkCellData>::iterator data2,
   vtkThreadLocal<vtkCellArray>::iterator data3,
   vtkThreadLocal<vtkCellArray>::iterator data4,
   vtkThreadLocal<vtkCellArray>::iterator data5,
   vtkThreadLocal<vtkCellArray>::iterator data6,
   vtkstd::vector<vtkIdType>::iterator offset1,
   vtkstd::vector<vtkIdType>::iterator offset2,
   vtkstd::vector<vtkIdType>::iterator offset3,
   vtkstd::vector<vtkIdType>::iterator offset4,
   vtkstd::vector<vtkIdType>::iterator offset5,
   vtkstd::vector<vtkIdType>::iterator offset6,
   vtkstd::vector<vtkIdType>::iterator offset7,
   vtkstd::vector<vtkIdType>::iterator offset8)
{
  if ( !this->MasterThreadPopulatedOutput )
    function->Execute( data1[0],
                       data2[0],
                       data3[0],
                       data4[0],
                       data5[0],
                       data6[0],
                       offset1[0],
                       offset2[0],
                       offset3[0],
                       offset4[0],
                       offset5[0],
                       offset6[0],
                       offset7[0],
                       offset8[0] );
}

void vtkRangeFunctor::ComputeMasterTID()
  {
  this->MasterThreadId = 0;
  }

void vtkTreeFunctor::ComputeMasterTID()
  {
  this->MasterThreadId = 0;
  }
