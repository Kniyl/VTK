#include "vtkSMP.h"

#include "vtkPoints.h"
#include "vtkFloatArray.h"
#include "vtkObjectFactory.h"

#include <time.h>
#include <sys/time.h>

//--------------------------------------------------------------------------------
vtkFunctor::vtkFunctor() { }

vtkFunctor::~vtkFunctor() { }

//--------------------------------------------------------------------------------
vtkFunctorInitialisable::vtkFunctorInitialisable() : vtkFunctor()
  {
  IsInitialized = false;
  }

vtkFunctorInitialisable::~vtkFunctorInitialisable() { }

bool vtkFunctorInitialisable::CheckAndSetInitialized() const
  {
  bool ret = IsInitialized;
  IsInitialized = true;
  return ret;
  }

//--------------------------------------------------------------------------------
vtkMergeable::vtkMergeable() { }

vtkMergeable::~vtkMergeable() { }

//--------------------------------------------------------------------------------
vtkMergeableInitialisable::vtkMergeableInitialisable() { }

vtkMergeableInitialisable::~vtkMergeableInitialisable() { }


namespace vtkSMP
{
  //--------------------------------------------------------------------------------
  void ForEach(vtkIdType first, vtkIdType last, const vtkFunctor& op)
    {
    InternalForEach( first, last, &op );
    }

  void ForEach(vtkIdType first, vtkIdType last, const vtkFunctorInitialisable& f )
    {
    if (!f.CheckAndSetInitialized())
      {
      struct timespec t0, t1;
      int ret_value = clock_gettime(CLOCK_REALTIME, &t0);
      InternalInit( &f );
      ret_value += clock_gettime(CLOCK_REALTIME, &t1);
      int s = t1.tv_sec - t0.tv_sec;
      int ns = t1.tv_nsec - t0.tv_nsec;
      if ( ns < 0 ) { s -= 1; ns += 1000000000; }
      if (ret_value) cout << "!";
      if (s)
        {
        cout << s;
        if ( ns < 100000000 ) cout << 0;
        if ( ns < 10000000 ) cout << 0;
        if ( ns < 1000000 ) cout << 0;
        if ( ns < 100000 ) cout << 0;
        if ( ns < 10000 ) cout << 0;
        if ( ns < 1000 ) cout << 0;
        if ( ns < 100 ) cout << 0;
        if ( ns < 10 ) cout << 0;
        }
      cout << ns << " ";
      }
    InternalForEach( first, last, &f );
    }

  void FillThreadsIDs( vtkstd::vector<vtkSMPThreadID>& result )
    {
    result.clear();
    InternalGetThreadsIDs( result );
    }

  void Merge( const vtkMergeable &f )
  {
    InternalMerge( &f );
  }

  void PreMerge( const vtkMergeableInitialisable &f )
  {
    InternalPreMerge( &f );
  }

  //--------------------------------------------------------------------------------
//  vtkStandardNewMacro(vtkThreadLocal);

//  vtkThreadLocal::vtkThreadLocal() : vtkObject() { }

//  vtkThreadLocal::~vtkThreadLocal()
//    {
//    for ( typename vtkstd::map<vtkSMPThreadID, vtkObject*>::iterator it = this->ThreadLocalStorage.begin();
//          it != this->ThreadLocalStorage.end(); ++it )
//      {
//      it->second->UnRegister(this);
//      it->second = 0;
//      }
//    this->ThreadLocalStorage.clear();
//    }

//  void vtkThreadLocal::PrintSelf(ostream &os, vtkIndent indent)
//    {
//    this->Superclass::PrintSelf( os, indent );
//    os << indent << "Class stored: ";
//    if (this->ThreadLocalStorage.size())
//      cout << this->ThreadLocalStorage[0]->GetClassName() << endl;
//    else
//      cout << "unknown" << endl;
//    os << indent << "Local storage: " << endl;
//    for ( vtkstd::map<vtkSMPThreadID, vtkObject*>::iterator it = this->ThreadLocalStorage.begin();
//          it != this->ThreadLocalStorage.end(); ++it )
//      {
//      os << indent.GetNextIndent() << "id " << it->first << ": (" << it->second << ")" << endl;
//      it->second->PrintSelf(os, indent.GetNextIndent().GetNextIndent());
//      }
//    }

//  void vtkThreadLocal::SetLocal ( vtkSMPThreadID tid, vtkObject* item )
//    {
//    if (this->ThreadLocalStorage[tid])
//      {
//      this->ThreadLocalStorage[tid]->UnRegister(this);
//      }

//    if (item)
//      {
//      item->Register(this);
//      }
//    this->ThreadLocalStorage[tid] = item;
//    }

//  vtkObject* vtkThreadLocal::GetLocal(vtkSMPThreadID tid)
//    {
//    return this->ThreadLocalStorage[tid];
//    }

}
