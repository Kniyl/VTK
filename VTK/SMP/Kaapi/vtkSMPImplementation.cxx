#include "vtkSMPImplementation.h"
#include <kaapic.h>

void smpInit(void)
{
  kaapic_init( 1 );
}
void smpFini(void)
{
  kaapic_finalize();
}

void functor( int32_t b, int32_t e, int32_t tid, const vtkFunctor* o, vtkIdType f )
{
  for ( int32_t k = b; k < e; ++k )
  {
    (*o)( f + k, tid );
  }
}

void InternalForEach ( vtkIdType first, vtkIdType last, const vtkFunctor* op )
{
  kaapic_foreach_attr_t attr;
  kaapic_foreach_attr_init(&attr);
  kaapic_foreach_attr_set_grains(&attr, 32, 32);
  kaapic_foreach( 0, last - first, &attr, 2, functor, op, first );
  kaapic_foreach_attr_destroy(&attr);
}

void InternalInit( const vtkFunctorInitialisable* f )
{
  vtkSMPThreadID numThreads = kaapic_get_concurrency();

  for (vtkSMPThreadID i = 0; i < numThreads; ++i)
  {
    f->init( i );
  }
}

void InternalGetThreadsIDs(vtkstd::vector<vtkSMPThreadID>& result)
{
  vtkSMPThreadID numThreads = kaapic_get_concurrency();

  for (vtkSMPThreadID i = 0; i < numThreads; ++i)
    result.push_back(i);
}
