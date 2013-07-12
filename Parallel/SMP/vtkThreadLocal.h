#ifndef __vtkThreadLocal_h__
#define __vtkThreadLocal_h__

#include "vtkParallelSMPModule.h" // For export macro
#include "vtkObject.h"
#include "vtkInstantiator.h"
#include <vector>
#include <typeinfo>

class vtkCellArray;
class vtkCellData;
class vtkParallelTree;
class vtkPointData;
class vtkPoints;
class vtkSMPMergePoints;
class vtkTask;

//======================================================================================
extern int vtkSMPInternalGetNumberOfThreads();

template<class T>
class VTKPARALLELSMP_EXPORT vtkThreadLocal : public vtkObject
{
  protected :
  vtkThreadLocal() :
      vtkObject(),
      ThreadLocalStorage(vtkSMPInternalGetNumberOfThreads(), NULL),
      SpecificName("NONE")
    {
    }
  ~vtkThreadLocal()
    {
    for ( iterator it = ThreadLocalStorage.begin();
          it != ThreadLocalStorage.end(); ++it )
      {
      if ( *it )
        (*it)->UnRegister( this );
      *it = 0;
      }
    ThreadLocalStorage.clear();
    }

  public:
  typedef typename std::vector<T*>::iterator iterator;

  vtkTypeMacro(vtkThreadLocal, vtkObject)
  static vtkThreadLocal<T>* New() { return new vtkThreadLocal<T>(); }

  void PrintSelf( ostream &os, vtkIndent indent )
    {
    this->Superclass::PrintSelf( os, indent );
    os << indent << "Class stored: " << typeid(T).name() << endl;
    os << indent << "Local storage: " << endl;
    size_t i = 0;
    for ( iterator it = ThreadLocalStorage.begin(); it != ThreadLocalStorage.end(); ++it, ++i )
      {
      os << indent.GetNextIndent() << "id " << i << ": (" << *it << ")" << endl;
      if ( *it ) (*it)->PrintSelf(os, indent.GetNextIndent().GetNextIndent());
      }
    }

  T* NewLocal ( T* specificImpl, int tid )
    {
    T* item = specificImpl->NewInstance();
    this->SetLocal(item,tid);

    if(item) item->Delete();

    return item;
    }

  T* NewLocal ( int tid )
    {
    if (this->SpecificName != "NONE")
      {
      return this->NewLocal(this->SpecificName.c_str(),tid);
      }

    T* item = T::New();
    this->SetLocal(item,tid);

    if(item) item->Delete();

    return item;
    }

  T* NewLocal ( const char* name, int tid )
    {
    vtkObject* obj = vtkInstantiator::CreateInstance(name);
    T* item = T::SafeDownCast(obj);

    if (item)
      {
      this->SetLocal(item,tid);
      item->Delete();
      return item;
      }

    this->SetLocal(0,tid);
    if(obj) obj->Delete();

    return 0;
    }

  iterator Begin( vtkIdType startItem = 0 )
    {
    iterator value = ThreadLocalStorage.begin();
    while ( startItem )
      {
      ++value;
      --startItem;
      }
    return value;
    }

  iterator End( )
    {
    return ThreadLocalStorage.end();
    }

  void SetSpecificClassName( const char* name )
    {
    this->SpecificName = std::string(name);
    }

  void SetLocal ( T* item, int tid )
    {
    if ( this->ThreadLocalStorage[tid] )
      {
      this->ThreadLocalStorage[tid]->UnRegister(this);
      }

    if ( item )
      {
      item->Register( this );
      }

    this->ThreadLocalStorage[tid] = item;
    }

  T*& GetLocal(int tid)
    {
    return this->ThreadLocalStorage[tid];
    }

  template<class Derived>
  Derived*& GetLocal(int tid)
    {
    return Derived::SafeDownCast(this->ThreadLocalStorage[tid]);
    }

  template<class Derived>
  void FillDerivedThreadLocal( vtkThreadLocal<Derived>* other )
    {
    T* elem;
    iterator src = ThreadLocalStorage.begin();
    for ( typename vtkThreadLocal<Derived>::iterator it = other->Begin();
          it != other->End(); ++it, ++src )
      {
      if ( (elem = *it) ) elem->UnRegister(other);
      Derived* d = (*it) = Derived::SafeDownCast(*src);
      if ( d ) d->Register(other);
      }
    other->SetSpecificClassName(this->SpecificName.c_str());
    }

 protected:
  std::vector<T*> ThreadLocalStorage;
  std::string SpecificName;
};

#endif //__vtkThreadLocal_h__
