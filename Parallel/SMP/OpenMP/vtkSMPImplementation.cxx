#include "vtkParallelOperators.h"
#include "vtkTreeFunctor.h"
#include "vtkTreeFunctorInitializable.h"
#include "vtkRangeFunctor.h"
#include "vtkRangeFunctorInitializable.h"
#include "vtkRange1D.h"
#include "vtkParallelTree.h"
#include "vtkTask.h"
#include "vtkMergeDataSets.h"

#include <omp.h>

void omp_traversal( vtkIdType index, int lvl, vtkIdType BranchingFactor, const vtkParallelTree* Tree, vtkTreeFunctor* func )
{
  int tid = omp_get_thread_num();
  vtkTreeFunctorInitializable* f = vtkTreeFunctorInitializable::SafeDownCast(func);
  if ( f && f->ShouldInitialize(tid) ) f->Init(tid);
  vtkLocalData* data = func->getLocal(tid);
  int cont = Tree->TraverseNode(index, lvl, func, data);
  data->Delete();
  if ( cont )
    {
    for ( vtkIdType i = index * BranchingFactor + 1, j = 0; j < BranchingFactor; ++i, ++j )
      {
#pragma omp task
      omp_traversal( i, lvl + 1, BranchingFactor, Tree, func );
      }
#pragma omp taskwait
    }
}

int vtkSMPInternalGetTid()
{
  return omp_get_thread_num();
}

int vtkSMPInternalGetNumberOfThreads()
{
  int numThreads = 1;
#pragma omp parallel shared(numThreads)
#pragma omp master
  {
  numThreads = omp_get_num_threads();
  }
  return numThreads;
}

void vtkParallelOperators::ForEach(vtkIdType first, vtkIdType last, const vtkRangeFunctor* op, int grain)
{
#pragma omp default(none)
#pragma omp parallel
  {
  vtkIdType chunkSize = (last - first) / omp_get_num_threads() + 1;
  int tid = omp_get_thread_num();
  vtkIdType f = first + chunkSize * tid;
  vtkIdType l = first + chunkSize * (tid + 1);
  if (l>last) l = last;
  vtkRange1D* range = vtkRange1D::New();
  range->Setup(f,l,tid);
  (*op)( range );
  range->Delete();
  }
}

void vtkParallelOperators::ForEach(vtkIdType first, vtkIdType last, const vtkRangeFunctorInitializable* op, int grain)
{
#pragma omp default(none)
#pragma omp parallel
  {
  vtkIdType chunkSize = (last - first) / omp_get_num_threads() + 1;
  int tid = omp_get_thread_num();
  if ( op->ShouldInitialize(tid) )
    op->Init(tid);
  vtkIdType f = first + chunkSize * tid;
  vtkIdType l = first + chunkSize * (tid + 1);
  if (l>last) l = last;
  vtkRange1D* range = vtkRange1D::New();
  range->Setup(f,l,tid);
  (*op)( range );
  range->Delete();
  }
}

void vtkMergeDataSets::Parallel(
    const vtkTask* function,
    vtkThreadLocal<vtkSMPMergePoints>::iterator data)
{
  int skipThreads = this->MasterThreadPopulatedOutput;
#pragma omp parallel shared(skipThreads, data)
  {
  int num_thread = omp_get_thread_num();
  if ( skipThreads <= num_thread )
    function->Execute( *(data+num_thread) );
  }
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
  int skipThreads = this->MasterThreadPopulatedOutput;
#pragma omp parallel shared(skipThreads, data1, data2, data3, data4, data5, data6, offset1, offset2, offset3, offset4, offset5, offset6, offset7, offset8 )
  {
  int num_thread = omp_get_thread_num();
  if ( skipThreads <= num_thread )
    function->Execute( *(data1+num_thread),
                       *(data2+num_thread),
                       *(data3+num_thread),
                       *(data4+num_thread),
                       *(data5+num_thread),
                       *(data6+num_thread),
                       *(offset1+num_thread),
                       *(offset2+num_thread),
                       *(offset3+num_thread),
                       *(offset4+num_thread),
                       *(offset5+num_thread),
                       *(offset6+num_thread),
                       *(offset7+num_thread),
                       *(offset8+num_thread) );
  }
}

void vtkParallelOperators::Traverse( const vtkParallelTree* Tree, vtkTreeFunctor* func )
{
  int lvl;
  vtkIdType bf;
  Tree->GetTreeSize(lvl,bf);
#pragma omp parallel
  {
#pragma omp single nowait
  {
  omp_traversal( 0, 0, bf, Tree, func);
  }
  }
}

void vtkRangeFunctor::ComputeMasterTID()
  {
  this->MasterThreadId = omp_get_thread_num();
  }

void vtkTreeFunctor::ComputeMasterTID()
  {
  this->MasterThreadId = omp_get_thread_num();
  }
