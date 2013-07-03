#include "vtkSMPImplementation.h"
#include "vtkParallelOperators.h"
#include "vtkFunctor.h"
#include "vtkFunctorInitializable.h"
#include "vtkParallelTree.h"
#include "vtkTask.h"
#include "vtkMergeDataSets.h"

int fillTid() {return -1;}

TBBInit::TBBInit() : init(), tids(fillTid)
  {
  atomique = -1;
  }

TBBInit::~TBBInit() {}

int TBBInit::getTID() const
  {
  int& tid = tids.local();
  if (tid == -1)
    {
    tid = ++atomique;
    }
  return tid;
  }

const TBBInit performInit;

class FuncCall
{
  const vtkFunctor* o;

public:
  void operator() ( const tbb::blocked_range<vtkIdType>& r ) const
    {
    vtkLocalData* data = o->getLocal(performInit.getTID());
    for ( vtkIdType k = r.begin(); k < r.end(); ++k )
      {
      (*o)( k, data );
      }
    data->Delete();
    }

  FuncCall ( const vtkFunctor* _o ) : o(_o) { }
  ~FuncCall () { }
};

class FuncCallInit
{
  const vtkFunctorInitializable* o;

public:
  void operator() ( const tbb::blocked_range<vtkIdType>& r ) const
    {
    int tid = performInit.getTID();
    if ( o->ShouldInitialize(tid) )
      {
      o->Init(tid);
      }
    vtkLocalData* data = o->getLocal(tid);
    for ( vtkIdType k = r.begin(); k < r.end(); ++k )
      {
      (*o)( k, data );
      }
    data->Delete();
    }

  FuncCallInit ( const vtkFunctorInitializable* _o ) : o(_o) { }
  ~FuncCallInit () { }
};

class TaskParallel : public tbb::task {
  public:
    const vtkTask* _task;
    vtkSMPMergePoints* _data;
    TaskParallel ( const vtkTask* t, vtkSMPMergePoints* d ) : _task(t), _data(d) {}
    tbb::task* execute()
      {
      _task->Execute(_data);
      return NULL;
      }
};

class TaskParallel2 : public tbb::task {
  public:
    const vtkTask* _task;
    vtkIdList* _data1;
    vtkCellData* _data2;
    vtkCellArray* _data3;
    vtkCellArray* _data4;
    vtkCellArray* _data5;
    vtkCellArray* _data6;
    vtkIdType _offset1, _offset2, _offset3, _offset4;
    vtkIdType _offset5, _offset6, _offset7, _offset8;

    TaskParallel2 ( const vtkTask* t, vtkIdList* d1, vtkCellData* d2,
        vtkCellArray* d3, vtkCellArray* d4, vtkCellArray* d5, vtkCellArray* d6,
        vtkIdType o1, vtkIdType o2, vtkIdType o3, vtkIdType o4,
        vtkIdType o5, vtkIdType o6, vtkIdType o7, vtkIdType o8 ) :
      _task(t), _data1(d1), _data2(d2), _data3(d3),
      _data4(d4), _data5(d5), _data6(d6),
      _offset1(o1), _offset2(o2), _offset3(o3), _offset4(o4),
      _offset5(o5), _offset6(o6), _offset7(o7), _offset8(o8) {}

    tbb::task* execute()
      {
      _task->Execute(
          _data1, _data2, _data3, _data4, _data5,
          _data6, _offset1, _offset2, _offset3,
          _offset4, _offset5, _offset6, _offset7, _offset8 );
      return NULL;
      }
};

class TaskTraverse : public tbb::task {
    const vtkParallelTree* Tree;
    vtkFunctor* Functor;
    const int level;
    const vtkIdType index, BranchingFactor;
    vtkLocalData* data;
  public:
    TaskTraverse(
        const vtkParallelTree* t, vtkFunctor* f,
        int l, vtkIdType i, vtkIdType b,
        vtkLocalData* d )
          : Tree(t), Functor(f),
            level(l), index(i),
            BranchingFactor(b),
            data(d)
          {
          data->Register(0);
          }
    ~TaskTraverse()
      {
      data->UnRegister(0);
      }

    tbb::task* execute()
      {
      if ( Tree->TraverseNode(index, level, Functor, data) )
        {
        int l = level + 1;
        this->set_ref_count(BranchingFactor + 1);
        TaskTraverse* t;
        for ( vtkIdType i = index * BranchingFactor + 1, j = 0;
              j < BranchingFactor; ++i, ++j )
          {
          t = new(this->allocate_child())
            TaskTraverse(Tree, Functor, l, i, BranchingFactor, data);
          tbb::task::spawn(*t);
          }
        this->wait_for_all();
        }
      return 0;
      }

    virtual void note_affinity ( tbb::task::affinity_id id )
      {
      int tid = id - 1;
      vtkFunctorInitializable* f =
        vtkFunctorInitializable::SafeDownCast(Functor);
      if ( f && f->ShouldInitialize(tid) ) f->Init(tid);
      data->UnRegister(0);
      data = Functor->getLocal(tid);
      data->Register(0);
      data->Delete();
      }
};

vtkIdType fillTLS() {return 0;}

int vtkSMPInternalGetNumberOfThreads()
{
  return tbb::task_scheduler_init::default_num_threads();
}

int vtkSMPInternalGetTid()
{
  return performInit.getTID();
}

void vtkParallelOperators::ForEach ( vtkIdType first, vtkIdType last, const vtkFunctor* op, int grain )
{
  vtkIdType n = last - first;
  if (!n) return;
  vtkIdType g = grain ? grain : sqrt(n);
  tbb::parallel_for( tbb::blocked_range<vtkIdType>( first, last, g ), FuncCall( op ) );
}

void vtkParallelOperators::ForEach ( vtkIdType first, vtkIdType last, const vtkFunctorInitializable* op, int grain )
{
  vtkIdType n = last - first;
  if (!n) return;
  vtkIdType g = grain ? grain : sqrt(n);
  tbb::parallel_for( tbb::blocked_range<vtkIdType>( first, last, g ), FuncCallInit( op ) );
}

void vtkMergeDataSets::Parallel(
    const vtkTask* function,
    vtkThreadLocal<vtkSMPMergePoints>::iterator data1)
  {
  int skipThreads = this->MasterThreadPopulatedOutput;
  for (vtkIdType tid = 0; tid < skipThreads; ++tid )
    {
    ++data1;
    }
  tbb::task_list list;
  for (vtkIdType tid = skipThreads; tid < tbb::task_scheduler_init::default_num_threads(); ++tid )
    {
    list.push_back(*new(tbb::task::allocate_root()) TaskParallel(function, *data1++));
    }
  tbb::task::spawn_root_and_wait(list);
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
  tbb::task_list list;
  for ( vtkIdType tid = skipThreads; tid < tbb::task_scheduler_init::default_num_threads(); ++tid )
    {
    list.push_back(*new(tbb::task::allocate_root()) TaskParallel2(
          function, *data1++, *data2++, *data3++, *data4++,
          *data5++, *data6++, *offset1++, *offset2++, *offset3++,
          *offset4++, *offset5++, *offset6++, *offset7++, *offset8++ ));
    }
  tbb::task::spawn_root_and_wait(list);
  }

/* */
void vtkParallelOperators::Traverse(const vtkParallelTree *Tree, vtkFunctor *func)
{
  int level;
  vtkIdType bf;
  Tree->GetTreeSize(level, bf);
  vtkLocalData* data = func->getLocal(performInit.getTID());
  TaskTraverse* t = new(tbb::task::allocate_root())
    TaskTraverse(Tree, func, 0, 0, bf, data);
  data->Delete();
  tbb::task::spawn_root_and_wait(*t);
}
/*/

struct ParallelTreeRange
{
  int root_level, cur_level, max;
  vtkIdType begin, end, *nodes_per_subtrees;

  ParallelTreeRange(int m, vtkIdType e, vtkIdType* n)
    : root_level(0), cur_level(0), max(m), begin(0), end(e),
    nodes_per_subtrees(n) {}

  bool empty() const {return begin == end; }
  bool is_divisible() const
  {
    return begin + nodes_per_subtrees[1] < end;
  }

  ParallelTreeRange(ParallelTreeRange& r, tbb::split)
  {
    begin = end = r.end;
    nodes_per_subtrees = r.nodes_per_subtrees;
    int steal_lvl = root_level + 1;
    while (r.cur_level >= steal_lvl && steal_lvl < r.max)
    {
      if (r.begin < r.end - r.nodes_per_subtrees[r.max - steal_lvl])
      {
        r.end = begin = end - r.nodes_per_subtrees[r.max-steal_lvl];
        cur_level = root_level = steal_lvl;
        max = r.max;
        break;
      }
      else
      {
        ++steal_lvl;
      }
    }
  }
};

vtkIdType convert_to_rowfirst(ParallelTreeRange& r, vtkIdType bf)
{
  if (!r.cur_level) return 0;
  int lvl = 1;
  vtkIdType index = 1;
  vtkIdType i = r.begin;
  while (lvl <= r.cur_level)
  {
    vtkIdType size = r.nodes_per_subtrees[r.max - lvl];
    while (i>size)
    {
      i -= size;
      ++index;
    }
    if (r.cur_level != lvl)
      index = index * bf + 1;
    ++lvl;
    --i;
  }
  return index;
}

struct Body
{
  const vtkParallelTree* Tree;
  vtkFunctor* func;
  vtkIdType branching_factor;

  Body(const vtkParallelTree* t, vtkFunctor* f, vtkIdType bf)
    : Tree(t), func(f), branching_factor(bf) {}

  void operator()(ParallelTreeRange& r) const
  {
    vtkIdType id = convert_to_rowfirst(r,branching_factor);
    int tid = performInit.getTID();
    vtkFunctorInitializable* init = vtkFunctorInitializable::SafeDownCast(func);
    if(init && init->ShouldInitialize(tid)) init->Init(tid);
    while(++r.begin<=r.end)
    {
      if (Tree->TraverseNode(id,r.cur_level,func,tid))
      {
        ++r.cur_level;
        id *= branching_factor;
      }
      else
      {
        if (r.max != r.cur_level)
        {
          r.begin += r.nodes_per_subtrees[r.max-r.cur_level]-1;
        }
        while ( !(id%branching_factor) && r.cur_level > r.root_level )
        {
          --r.cur_level;
          id = (id-1)/branching_factor;
        }
      }
      ++id;
    }
  }
};

void vtkParallelOperators::Traverse(const vtkParallelTree *Tree, vtkFunctor *func)
{
  int max_level;
  vtkIdType branching_factor;
  Tree->GetTreeSize(max_level,branching_factor);
  Body b(Tree, func, branching_factor);

  vtkIdType *npst = new vtkIdType[max_level + 1];
  vtkIdType value = npst[0] = 1;
  for (int i = 1; i <= max_level; ++i)
  {
    value *= branching_factor;
    npst[i] = ++value;
  }
  ParallelTreeRange r(max_level,value,npst);

  tbb::parallel_for(r,b);
  delete npst;
}
/* */

void vtkFunctor::ComputeMasterTID()
  {
  this->MasterThreadId = performInit.getTID();
  }
