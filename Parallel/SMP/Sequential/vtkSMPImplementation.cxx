#include "vtkParallelOperators.h"
#include "vtkFunctor.h"
#include "vtkFunctorInitializable.h"
#include "vtkParallelTree.h"
#include "vtkTask.h"
#include "vtkMergeDataSets.h"

void sequential_traverse( vtkIdType index, int lvl, vtkIdType BranchingFactor, const vtkParallelTree* Tree, vtkFunctor* func, vtkLocalData* data )
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

void vtkParallelOperators::ForEach( vtkIdType first, vtkIdType last, const vtkFunctor* op, int grain )
{
  vtkLocalData* data = op->getLocal(0);
  for ( ; first < last; ++first )
    (*op)( first, data );
  data->Delete();
}

void vtkParallelOperators::ForEach( vtkIdType first, vtkIdType last, const vtkFunctorInitializable* f, int grain )
{
  if ( f->ShouldInitialize(0) )
    f->Init(0);
  vtkLocalData* data = f->getLocal(0);
  for ( ; first < last; ++first )
    (*f)( first, data );
  data->Delete();
}

void vtkParallelOperators::Traverse( const vtkParallelTree* Tree, vtkFunctor* func )
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

void vtkFunctor::ComputeMasterTID()
  {
  this->MasterThreadId = 0;
  }
