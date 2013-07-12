#include "vtkSMPImplementation.h"
#include "vtkParallelOperators.h"
#include "vtkTreeFunctor.h"
#include "vtkTreeFunctorInitializable.h"
#include "vtkRangeFunctor.h"
#include "vtkRangeFunctorInitializable.h"
#include "vtkRange1D.h"
#include "vtkParallelTree.h"
#include "vtkTask.h"
#include "vtkMergeDataSets.h"
#include <kaapic.h>
#include <cmath>

KaapiInit::KaapiInit() // : com(ka::System::join_community()) { }
  {
  kaapic_init(KAAPIC_START_ONLY_MAIN);
  }

KaapiInit::~KaapiInit()
  {
//  com.leave();
//  ka::System::terminate();
  kaapic_finalize();
  }

const KaapiInit performInit;

struct TaskParallel2 : public ka::Task<15>::Signature<
        const vtkTask*, vtkIdList*, vtkCellData*, vtkCellArray*, vtkCellArray*,
        vtkCellArray*, vtkCellArray*, vtkIdType, vtkIdType, vtkIdType,
        vtkIdType, vtkIdType, vtkIdType, vtkIdType, vtkIdType> {};

template<> struct TaskBodyCPU<TaskParallel2> {
  void operator() ( const vtkTask* function,
                    vtkIdList* data1, vtkCellData* data2,
                    vtkCellArray* data3, vtkCellArray* data4,
                    vtkCellArray* data5, vtkCellArray* data6,
                    vtkIdType offset1, vtkIdType offset2,
                    vtkIdType offset3, vtkIdType offset4,
                    vtkIdType offset5, vtkIdType offset6,
                    vtkIdType offset7, vtkIdType offset8 )
    {
    function->Execute(
        data1, data2, data3, data4, data5,
        data6, offset1, offset2, offset3,
        offset4, offset5, offset6, offset7, offset8 );
    }
};

struct TaskParallel : public ka::Task<2>::Signature<const vtkTask*, vtkSMPMergePoints*> {};

template<> struct TaskBodyCPU<TaskParallel> {
  void operator() ( const vtkTask* function, vtkSMPMergePoints* data )
    {
    function->Execute(data);
    }
};

typedef struct tree_work
{
  kaapi_workqueue_t wq;
  const vtkParallelTree* Tree;
  vtkTreeFunctor* op;
  int max_level;
  int cur_level;
  vtkIdType branching_factor;
  vtkIdType* nodes_per_subtrees;
  int root_lvl;
} work_t;

vtkIdType convert_to_row_first( kaapi_workqueue_index_t i, int level, work_t* work )
  {
  if ( !level ) return 0;
  int lvl = 1;
  vtkIdType index = 1;
  while ( lvl <= level )
    {
    vtkIdType size = work->nodes_per_subtrees[work->max_level - lvl];
    while ( i > size )
      {
      i -= size;
      ++index;
      }
    if ( level != lvl )
      index = index * work->branching_factor + 1;
    ++lvl;
    --i;
    }
  return index;
  }

static void thief_entrypoint( void* args, kaapi_thread_t* thread )
  {
  work_t* const work = (work_t*)(args);

  int tid = kaapi_get_self_kid();
  vtkTreeFunctorInitializable* iop = vtkTreeFunctorInitializable::SafeDownCast( work->op );
  if ( iop && iop->ShouldInitialize(tid) )
    iop->Init(tid);
  vtkLocalData* data = work->op->getLocal(tid);

  kaapi_workqueue_index_t i, nil;

  vtkIdType id = convert_to_row_first(kaapi_workqueue_range_begin(&work->wq), work->cur_level, work);
  while ( !kaapi_workqueue_pop(&work->wq, &i, &nil, 1) )
    {
    if ( work->Tree->TraverseNode( id, work->cur_level, work->op, data ) )
      {
      ++work->cur_level;
      id *= work->branching_factor;
      }
    else
      {
      if ( work->max_level != work->cur_level )
        {
        kaapi_workqueue_pop(&work->wq, &i, &nil, work->nodes_per_subtrees[work->max_level - work->cur_level] - 1);
        }
      while ( !(id % work->branching_factor) && work->cur_level > work->root_lvl )
        {
        --work->cur_level;
        id = (id - 1) / work->branching_factor;
        }
      }
    ++id;
    }

  data->Delete();

  }

static int splitter(
    struct kaapi_task_t* victim_task,
    void*  args,
    struct kaapi_listrequest_t* lr,
    struct kaapi_listrequest_iterator_t* lri
    )
  {
  work_t* const work = (work_t*)(args);
  kaapi_workqueue_index_t i, j;
  int steal_lvl = work->root_lvl + 1;

  kaapi_request_t* req = kaapi_api_listrequest_iterator_get(lr, lri);
  while ( req && work->cur_level >= steal_lvl && steal_lvl < work->max_level )
    {
    if ( !(kaapi_workqueue_steal(&work->wq, &i, &j, work->nodes_per_subtrees[work->max_level - steal_lvl]) ) )
      {
      work_t* const tw = (work_t*)kaapi_request_pushdata(req, sizeof(work_t));
      tw->Tree = work->Tree;
      tw->op = work->op;
      tw->branching_factor = work->branching_factor;
      tw->max_level = work->max_level;
      tw->nodes_per_subtrees = work->nodes_per_subtrees;
      tw->root_lvl = steal_lvl;
      tw->cur_level = steal_lvl;
      kaapi_workqueue_init_with_kproc( &tw->wq, i, j, req->ident );
      kaapi_task_init( kaapi_request_toptask(req), thief_entrypoint, tw);
      kaapi_request_pushtask_adaptive( req, victim_task, splitter, 0 );
      kaapi_request_committask(req);

      req = kaapi_api_listrequest_iterator_next(lr, lri);
      }
    else
      ++steal_lvl;
    }

  return 0;
  }

inline void doFor( int32_t b, int32_t e, int32_t tid, const vtkRangeFunctor* o )
  {
  vtkRange1D* range = vtkRange1D::New();
  range->Setup(b,e,tid);
  (*o)(range);
  range->Delete();
  }
inline void doForInit( int32_t b, int32_t e, int32_t tid, const vtkRangeFunctorInitializable* o )
  {
  if (o->ShouldInitialize(tid))
    o->Init(tid);
  doFor(b,e,tid,o);
  }

//--------------------------------------------------------------------------------
int vtkSMPInternalGetNumberOfThreads()
{
  return kaapi_getconcurrency();
}

int vtkSMPInternalGetTid()
{
  return kaapi_get_self_kid();
}

void vtkParallelOperators::ForEach ( vtkIdType first, vtkIdType last, const vtkRangeFunctor* op, int grain )
{
  vtkIdType n = last - first;
  int g = grain ? grain : sqrt(n);
  kaapic_begin_parallel(KAAPIC_FLAG_DEFAULT);
  kaapic_foreach_attr_t attr;
  kaapic_foreach_attr_init(&attr);
  kaapic_foreach_attr_set_grains(&attr, g, g);
  kaapic_foreach( first, last, &attr, 1, doFor, op );
  kaapic_end_parallel(KAAPIC_FLAG_DEFAULT);
  kaapic_foreach_attr_destroy(&attr);
}

void vtkParallelOperators::ForEach ( vtkIdType first, vtkIdType last, const vtkRangeFunctorInitializable* op, int grain )
{
  vtkIdType n = last - first;
  int g = grain ? grain : sqrt(n);
  kaapic_begin_parallel(KAAPIC_FLAG_DEFAULT);
  kaapic_foreach_attr_t attr;
  kaapic_foreach_attr_init(&attr);
  kaapic_foreach_attr_set_grains(&attr, g, g);
  kaapic_foreach( first, last, &attr, 1, doForInit, op );
  kaapic_end_parallel(KAAPIC_FLAG_DEFAULT);
  kaapic_foreach_attr_destroy(&attr);
}

void vtkMergeDataSets::Parallel(
    const vtkTask* function,
    vtkThreadLocal<vtkSMPMergePoints>::iterator data1)
{
  int skipThreads = this->MasterThreadPopulatedOutput;
  for ( vtkIdType tid = 0; tid < skipThreads; ++tid )
    {
    ++data1;
    }
//    kaapi_begin_parallel( KAAPI_SCHEDFLAG_DEFAULT );
  kaapic_begin_parallel(KAAPIC_FLAG_DEFAULT);
  for ( vtkIdType tid = skipThreads; tid < kaapi_getconcurrency(); ++tid )
    {
    ka::Spawn<TaskParallel>()( function, *data1++ );
    }
  kaapic_end_parallel(KAAPIC_FLAG_DEFAULT);
//    kaapi_end_parallel( KAAPI_SCHEDFLAG_DEFAULT );
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
  for ( vtkIdType tid = 0; tid < skipThreads; ++tid )
    {
    ++data1; ++data2; ++data3; ++data4; ++data5; ++data6;
    ++offset1; ++offset2; ++offset3; ++offset4; ++offset5; ++offset6; ++offset7; ++offset8;
    }
//    kaapi_begin_parallel( KAAPI_SCHEDFLAG_DEFAULT );
  kaapic_begin_parallel(KAAPIC_FLAG_DEFAULT);
  for ( vtkIdType tid = skipThreads; tid < kaapi_getconcurrency(); ++tid )
    {
    ka::Spawn<TaskParallel2>()(
        function, *data1++, *data2++, *data3++, *data4++, *data5++, *data6++,
        *offset1++, *offset2++, *offset3++, *offset4++,
        *offset5++, *offset6++, *offset7++, *offset8++ );
    }
  kaapic_end_parallel(KAAPIC_FLAG_DEFAULT);
//    kaapi_end_parallel( KAAPI_SCHEDFLAG_DEFAULT );
}

void vtkParallelOperators::Traverse( const vtkParallelTree *Tree, vtkTreeFunctor* func )
{
  work_t work;
  kaapi_workqueue_index_t i, nil;

  kaapic_begin_parallel(KAAPIC_FLAG_DEFAULT);

  /* get the self thread */
  kaapi_thread_t* thread = kaapi_self_thread();

  /* initialize work */
  work.Tree = Tree;
  work.op = func;
  Tree->GetTreeSize( work.max_level, work.branching_factor );
  work.nodes_per_subtrees = new vtkIdType[work.max_level + 1];
  vtkIdType value = work.nodes_per_subtrees[work.root_lvl = 0] = 1;
  for ( vtkIdType i = 1; i <= work.max_level; ++i )
    {
    value *= work.branching_factor;
    work.nodes_per_subtrees[i] = ++value;
    }
  kaapi_workqueue_init(&work.wq, 0, (kaapi_workqueue_index_t)value);

  /* push an adaptive task */
  void* sc = kaapi_task_begin_adaptive(
                                       thread,
                                       KAAPI_SC_CONCURRENT | KAAPI_SC_NOPREEMPTION,
                                       splitter,
                                       &work     /* arg for splitter = work to split */
                                       );

  work.cur_level = work.root_lvl;
  vtkLocalData* data = func->getLocal(kaapi_get_self_kid());
  vtkIdType id = 0;
  while ( !kaapi_workqueue_pop(&work.wq, &i, &nil, 1) )
    {
    if ( Tree->TraverseNode( id, work.cur_level, func, data ) )
      {
      ++work.cur_level;
      id *= work.branching_factor;
      }
    else
      {
      if ( work.max_level != work.cur_level )
        {
        kaapi_workqueue_pop(&work.wq, &i, &nil, work.nodes_per_subtrees[work.max_level - work.cur_level] - 1);
        }
      while ( !(id % work.branching_factor) && work.cur_level > work.root_lvl )
        {
        --work.cur_level;
        id = (id - 1) / work.branching_factor;
        }
      }
    ++id;
    }
  data->Delete();

  kaapi_task_end_adaptive(thread, sc);

  /* wait for thieves */
  kaapi_sched_sync( );

  delete [] work.nodes_per_subtrees;

  kaapic_end_parallel(KAAPIC_FLAG_DEFAULT);

}

void vtkTreeFunctor::ComputeMasterTID()
  {
  this->MasterThreadId = kaapi_get_self_kid();
  }

void vtkRangeFunctor::ComputeMasterTID()
  {
  this->MasterThreadId = kaapi_get_self_kid();
  }
