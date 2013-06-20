#ifndef __vtkSMP_h__
#define __vtkSMP_h__

#include "vtkSMPImplementation.h"
#include "vtkInstantiator.h"
#include "vtkObject.h"
#include "vtkPoints.h"
#include "vtkPointData.h"
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkIdList.h"
#include "vtkSMPMergePoints.h"
#include <vector>
#include <typeinfo>

namespace vtkSMP
{
  class vtkIdTypeThreadLocal;
  template <class T> class vtkThreadLocal;
  int InternalGetNumberOfThreads();
  int InternalGetTid();
}

template <typename... Params>
class VTK_SMP_EXPORT vtkFunctor : public vtkObject
{
  vtkFunctor ( const vtkFunctor& );
  void operator =( const vtkFunctor& );

public:
  vtkTypeMacro(vtkFunctor,vtkObject);
  void PrintSelf(ostream& os, vtkIndent indent)
    {
    this->Superclass::PrintSelf(os,indent);
    }

  virtual void operator () ( Params... params ) const = 0;

protected:
  vtkFunctor(){}
  ~vtkFunctor(){}
};

template <typename... Params>
class VTK_SMP_EXPORT vtkFancyFunctor : public vtkFunctor<>
{
  vtkFancyFunctor(const vtkFancyFunctor&);
  void operator =(const vtkFancyFunctor&);
public:
  vtkTypeMacro(vtkFancyFunctor,vtkFunctor<>);
  void PrintSelf(ostream &os, vtkIndent indent)
    {
    this->Superclass::PrintSelf(os,indent);
    }
  virtual void ThreadedMoveBasePointer(Params...) const = 0;
protected:
  vtkFancyFunctor(){}
  ~vtkFancyFunctor(){}
};

template <typename... Params>
class VTK_SMP_EXPORT vtkFunctorInitialisable : public vtkFunctor<Params...>
{
  vtkFunctorInitialisable ( const vtkFunctorInitialisable& );
  void operator =( const vtkFunctorInitialisable& );

public:
  vtkTypeMacro(vtkFunctorInitialisable,vtkFunctor<Params...>);
  void PrintSelf(ostream& os, vtkIndent indent)
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

  virtual void Init ( ) const = 0;
  bool ShouldInitialize( ) const
    {
    return !IsInitialized[vtkSMP::InternalGetTid()];
    }

protected:
  void Initialized( ) const
    {
    IsInitialized[vtkSMP::InternalGetTid()] = 1;
    }
  mutable vtkstd::vector<vtkIdType> IsInitialized;

  vtkFunctorInitialisable() : vtkFunctor<Params...>(),
    IsInitialized(vtkSMP::InternalGetNumberOfThreads(),0) {}
  ~vtkFunctorInitialisable(){ IsInitialized.clear(); }
};

template <class T> struct TaskHelper { typedef T* type; };
template<> struct TaskHelper<vtkIdType> { typedef vtkIdType type; };

template <class... T>
class VTK_SMP_EXPORT vtkTask : public vtkObjectBase
{
  vtkTask(const vtkTask&);
  void operator =(const vtkTask&);

public:
  vtkTypeMacro(vtkTask<T...>, vtkObjectBase);
  void PrintSelf(ostream &os, vtkIndent indent)
    {
    this->Superclass::PrintSelf(os,indent);
    }

  virtual void Execute(typename TaskHelper<T>::type... args) const = 0;

protected:
  vtkTask(){}
  ~vtkTask(){}
};

template <>
class VTK_SMP_EXPORT vtkTask<> : public vtkObjectBase
{
  vtkTask(const vtkTask&);
  void operator =(const vtkTask&);

public:
  vtkTypeMacro(vtkTask<>, vtkObjectBase);
  void PrintSelf(ostream &os, vtkIndent indent)
    {
    this->Superclass::PrintSelf(os,indent);
    }

  virtual void Execute() = 0;

protected:
  vtkTask(){}
  ~vtkTask(){}
};

class VTK_SMP_EXPORT vtkParallelTree
{
public:
  virtual int TraverseNode( vtkIdType id, int lvl, vtkFunctor<vtkIdType>* function ) const = 0;
  virtual void GetTreeSize ( int& max_level, vtkIdType& branching_factor ) const = 0;
};

template <class T> struct IteratorHelper { typedef typename vtkSMP::vtkThreadLocal<T>::iterator type; };
template<> struct IteratorHelper<vtkIdType> { typedef vtkstd::vector<vtkIdType>::iterator type; };

#include "vtkSMPImplementation.txx"

namespace vtkSMP
{
  void VTK_SMP_EXPORT ThreadedPrint(const char* message);

  template<bool v, class T>
  class vtkThreadLocalImpl : public vtkObject
    {
    protected :
      vtkThreadLocalImpl() : vtkObject(), ThreadLocalStorage(InternalGetNumberOfThreads(), NULL), Specific("NONE") {}
      ~vtkThreadLocalImpl()
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
      typedef typename vtkstd::vector<T*>::iterator iterator;

      vtkTypeMacro(vtkThreadLocalImpl, vtkObject)

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

      T* NewLocal ( T* specificImpl )
        {
        T* item = specificImpl->NewInstance();
        this->SetLocal( item );
        if(item) item->Delete();

        return item;
        }

      T* NewLocal ( )
        {
        if ( this->Specific != "NONE" )
          {
          return this->NewLocal(this->Specific.c_str());
          }
        T* item = static_cast<T*>(T::New());
        this->SetLocal( item );
        if(item) item->Delete();

        return item;
        }

      T* NewLocal ( const char* type )
        {
        vtkObject* obj = vtkInstantiator::CreateInstance(type);
        T* item = T::SafeDownCast(obj);
        if (item)
          {
          this->SetLocal(item);
          item->Delete();
          return item;
          }
        this->SetLocal(0);
        if (obj) obj->Delete();
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

      void SetSpecificClassName( const char* type )
        {
        this->Specific = vtkstd::string(type);
        }

      void SetLocal ( T* item )
        {
        int tid = InternalGetTid();
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

      T* GetLocal()
        {
        return this->ThreadLocalStorage[InternalGetTid()];
        }

      template<class Derived>
      Derived* GetLocal()
        {
        return Derived::SafeDownCast(this->ThreadLocalStorage[InternalGetTid()]);
        }

      template<class Derived>
      void FillDerivedThreadLocal( vtkThreadLocalImpl<v, Derived>* other )
        {
        T* elem;
        iterator src = ThreadLocalStorage.begin();
        for ( typename vtkSMP::vtkThreadLocalImpl<v, Derived>::iterator it = other->Begin();
              it != other->End(); ++it, ++src )
          {
          if ( (elem = *it) ) elem->UnRegister(other);
          Derived* d = (*it) = Derived::SafeDownCast(*src);
          if ( d ) d->Register(other);
          }
        other->SetSpecificClassName(this->Specific.c_str());
        }

    protected:
      vtkstd::vector<T*> ThreadLocalStorage;
      vtkstd::string Specific;
    };

  template<class T>
  class vtkThreadLocalImpl<false, T> : public vtkObject
    {
    protected :
      vtkThreadLocalImpl() : vtkObject(), ThreadLocalStorage(InternalGetNumberOfThreads(), NULL) {}
      ~vtkThreadLocalImpl()
        {
        ThreadLocalStorage.clear();
        }

    public:
      typedef typename vtkstd::vector<T*>::iterator iterator;

      vtkTypeMacro(vtkThreadLocalImpl, vtkObject)

      void PrintSelf( ostream &os, vtkIndent indent )
        {
        this->Superclass::PrintSelf( os, indent );
        os << indent << "Class stored: " << typeid(T).name() << endl;
        os << indent << "Local storage: " << endl;
        size_t i = 0;
        for ( iterator it = ThreadLocalStorage.begin(); it != ThreadLocalStorage.end(); ++it, ++i )
          {
          os << indent.GetNextIndent() << "id " << i << ": (" << *it << ") " << *(*it) << endl;
          }
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

      void SetLocal ( T* item )
        {
        int tid = InternalGetTid();
        this->ThreadLocalStorage[tid] = item;
        }

      T* GetLocal()
        {
        return this->ThreadLocalStorage[InternalGetTid()];
        }

      template<class Derived>
      Derived* GetLocal()
        {
        return static_cast<Derived*>(this->ThreadLocalStorage[InternalGetTid()]);
        }

      template<class Derived>
      void FillDerivedThreadLocal( vtkThreadLocalImpl<false, Derived>* other )
        {
        iterator src = ThreadLocalStorage.begin();
        for ( typename vtkSMP::vtkThreadLocalImpl<false, Derived>::iterator it = other->Begin();
              it != other->End(); ++it, ++src )
          {
          (*it) = static_cast<Derived*>(*src);
          }
        }

    protected:
      vtkstd::vector<T*> ThreadLocalStorage;
    };

  template <class T>
  class VTK_SMP_EXPORT vtkThreadLocal : public vtkThreadLocalImpl<std::is_base_of<vtkObject,T>::value, T>
  {
    vtkThreadLocal(const vtkThreadLocal&);
    void operator =(const vtkThreadLocal&);
  public:
    static vtkThreadLocal<T>* New() { return new vtkThreadLocal<T>(); }
#define MY_COMMA ,
    vtkTypeMacro(vtkThreadLocal, vtkThreadLocalImpl<std::is_base_of<vtkObject MY_COMMA T>::value MY_COMMA T>);
#undef MY_COMMA
    void PrintSelf(ostream &os, vtkIndent indent)
      {
      this->Superclass::PrintSelf(os,indent);
      }
  protected:
    vtkThreadLocal(){}
    ~vtkThreadLocal(){}
  };

  // ForEach template : parallel loop over an iterator
  void VTK_SMP_EXPORT ForEach( vtkIdType first, vtkIdType last, const vtkFunctor<vtkIdType>* op, int grain = 0 );

  void VTK_SMP_EXPORT ForEach( vtkIdType first, vtkIdType last, const vtkFunctorInitialisable<vtkIdType>* f, int grain = 0 );

  void VTK_SMP_EXPORT ForEach( vtkIdType first0, vtkIdType last0, vtkIdType first1, vtkIdType last1, const vtkFunctor<vtkIdType,vtkIdType>* op, int grain = 0 );

  void VTK_SMP_EXPORT ForEach( vtkIdType first0, vtkIdType last0, vtkIdType first1, vtkIdType last1, const vtkFunctorInitialisable<vtkIdType,vtkIdType>* f, int grain = 0 );

  void VTK_SMP_EXPORT ForEach( vtkIdType first0, vtkIdType last0, vtkIdType first1, vtkIdType last1, vtkIdType first2, vtkIdType last2, const vtkFunctor<vtkIdType,vtkIdType,vtkIdType>* op, int grain = 0 );

  void VTK_SMP_EXPORT ForEach( vtkIdType first0, vtkIdType last0, vtkIdType first1, vtkIdType last1, vtkIdType first2, vtkIdType last2, const vtkFunctorInitialisable<vtkIdType,vtkIdType,vtkIdType>* f, int grain = 0 );

  void VTK_SMP_EXPORT StaticForEach( vtkIdType first, vtkIdType last, const vtkFunctor<vtkIdType>* op );

  void VTK_SMP_EXPORT StaticForEach( vtkIdType first, vtkIdType last, const vtkFunctorInitialisable<vtkIdType>* op );

  void VTK_SMP_EXPORT ForEach( vtkIdType first0, vtkIdType last0, vtkIdType first1, vtkIdType last1, const vtkFancyFunctor<vtkIdType,vtkIdType>* op, int grain = 0 );

  void VTK_SMP_EXPORT ForEach( vtkIdType first0, vtkIdType last0, vtkIdType first1, vtkIdType last1, vtkIdType first2, vtkIdType last2, const vtkFancyFunctor<vtkIdType,vtkIdType,vtkIdType>* op, int grain = 0 );

  template <class... T>
  void VTK_SMP_EXPORT Parallel( const vtkTask<T...>* function, vtkIdType skipThreads, typename IteratorHelper<T>::type... data )
    {
    if (skipThreads)
      ParallelSkip<T...>( function, 0, skipThreads, data... );
    else
      ParallelDo<T...>( function, 0, InternalGetNumberOfThreads(), data... );
    }

  void VTK_SMP_EXPORT Traverse( const vtkParallelTree* Tree, vtkFunctor<vtkIdType>* func );

  void VTK_SMP_EXPORT MergePoints( vtkPoints* outPoints, vtkThreadLocal<vtkPoints>* inPoints, const double bounds[6],
                                   vtkPointData* outPtsData, vtkThreadLocal<vtkPointData>* inPtsData,
                                   vtkCellArray* outVerts, vtkThreadLocal<vtkCellArray>* inVerts,
                                   vtkCellArray* outLines, vtkThreadLocal<vtkCellArray>* inLines,
                                   vtkCellArray* outPolys, vtkThreadLocal<vtkCellArray>* inPolys,
                                   vtkCellArray* outStrips, vtkThreadLocal<vtkCellArray>* inStrips,
                                   vtkCellData* outCellsData, vtkThreadLocal<vtkCellData>* inCellsData, int SkipThreads );

  void VTK_SMP_EXPORT MergePoints( vtkSMPMergePoints* outPoints, vtkThreadLocal<vtkSMPMergePoints>* inPoints,
                                   vtkPointData* outPtsData, vtkThreadLocal<vtkPointData>* inPtsData,
                                   vtkCellArray* outVerts, vtkThreadLocal<vtkCellArray>* inVerts,
                                   vtkCellArray* outLines, vtkThreadLocal<vtkCellArray>* inLines,
                                   vtkCellArray* outPolys, vtkThreadLocal<vtkCellArray>* inPolys,
                                   vtkCellArray* outStrips, vtkThreadLocal<vtkCellArray>* inStrips,
                                   vtkCellData* outCellsData, vtkThreadLocal<vtkCellData>* inCellsData, int SkipThreads );

  class VTK_SMP_EXPORT vtkSpawnTasks : public vtkObject {
    protected:
      vtkSpawnTasks();
      ~vtkSpawnTasks();

      void InternalSpawn( vtkTask<>* );

    public:
      vtkTypeMacro(vtkSpawnTasks,vtkObject);
      static vtkSpawnTasks* New();
      void PrintSelf(ostream &os, vtkIndent indent);

    private:
      vtkSpawnTasks( const vtkSpawnTasks& );
      void operator =( const vtkSpawnTasks& );
  };
}

#endif //__vtkSMP_h__
