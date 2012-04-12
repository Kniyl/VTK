#ifndef __vtkSMP_h__
#define __vtkSMP_h__

#include "vtkSMPImplementation.h"
#include "vtkObject.h"
#include "vtkCommand.h"

#include <vector>
#include <typeinfo>  //for 'typeid'

class vtkPoints;
class vtkCommand;
class vtkPointData;
class vtkCellArray;
class vtkCellData;
class vtkSMPMergePoints;

class VTK_SMP_EXPORT vtkFunctor : public vtkObject
{
  vtkFunctor ( const vtkFunctor& );
  void operator =( const vtkFunctor& );

public:
  vtkTypeMacro(vtkFunctor,vtkObject);
  void PrintSelf(ostream& os, vtkIndent indent);

  virtual void operator () ( vtkIdType, vtkSMPThreadID ) const = 0;

protected:
  vtkFunctor();
  ~vtkFunctor();
};

class VTK_SMP_EXPORT vtkFunctorInitialisable : public vtkFunctor
{
  vtkFunctorInitialisable ( const vtkFunctorInitialisable& );
  void operator =( const vtkFunctorInitialisable& );

public:
  vtkTypeMacro(vtkFunctorInitialisable,vtkFunctor);
  void PrintSelf(ostream& os, vtkIndent indent);

  virtual void Init ( vtkSMPThreadID ) const = 0;
  bool ShouldInitialize( vtkSMPThreadID tid ) const;

protected:
  mutable vtkIdType* IsInitialized;

  vtkFunctorInitialisable();
  ~vtkFunctorInitialisable();
};

class VTK_SMP_EXPORT vtkTask : public vtkObjectBase
{
  vtkTask(const vtkTask&);
  void operator =(const vtkTask&);

public:
  vtkTypeMacro(vtkTask, vtkObjectBase);
  void PrintSelf(ostream &os, vtkIndent indent);

  virtual void Execute( vtkSMPThreadID tid, const vtkObject *data ) const = 0;

protected:
  vtkTask();
  ~vtkTask();
};

class VTK_SMP_EXPORT vtkTaskSplitable : public vtkTask
{
  vtkTaskSplitable ( const vtkTaskSplitable& );
  void operator =( const vtkTaskSplitable& );

public:
  vtkTypeMacro(vtkTaskSplitable, vtkTask);
  void PrintSelf(ostream &os, vtkIndent indent);

  virtual void SplitTask() const = 0;

protected:
  vtkTaskSplitable();
  ~vtkTaskSplitable();
};

namespace vtkSMP
{
  vtkSMPThreadID VTK_SMP_EXPORT GetNumberOfThreads( );


  template<class T>
  class VTK_SMP_EXPORT vtkThreadLocal : public vtkObject
    {
    protected :
      vtkThreadLocal() : vtkObject(), ThreadLocalStorage(GetNumberOfThreads()) { }
      ~vtkThreadLocal()
        {
        for ( typename vtkstd::vector<T*>::iterator it = ThreadLocalStorage.begin();
              it != ThreadLocalStorage.end(); ++it )
          {
          (*it)->UnRegister( this );
          *it = 0;
          }
        ThreadLocalStorage.clear();
        }

    public:
      vtkTypeMacro(vtkThreadLocal, vtkObject)
      static vtkThreadLocal<T>* New()
        {
        return new vtkThreadLocal<T>();
        }

      void PrintSelf( ostream &os, vtkIndent indent )
        {
        this->Superclass::PrintSelf( os, indent );
        os << indent << "Class stored: " << typeid(T).name() << endl;
        os << indent << "Local storage: " << endl;
        for ( size_t i = 0; i < this->ThreadLocalStorage.size(); ++i )
          {
          os << indent.GetNextIndent() << "id " << i << ": (" << ThreadLocalStorage[i] << ")" << endl;
          ThreadLocalStorage[i]->PrintSelf(os, indent.GetNextIndent().GetNextIndent());
          }
        }

      T* NewLocal ( vtkSMPThreadID tid, T* specificImpl )
        {
        if (this->ThreadLocalStorage[tid])
          {
          this->ThreadLocalStorage[tid]->UnRegister(this);
          }

        T* item = specificImpl->NewInstance();
        if (item)
          {
          item->Register(this);
          item->Delete();
          }
        this->ThreadLocalStorage[tid] = item;

        return item;
        }

      T* NewLocal ( vtkSMPThreadID tid )
        {
        if (this->ThreadLocalStorage[tid])
          {
          this->ThreadLocalStorage[tid]->UnRegister(this);
          }

        T* item = T::New();
        if (item)
          {
          item->Register(this);
          item->Delete();
          }
        this->ThreadLocalStorage[tid] = item;

        return item;
        }

      void SetLocal ( vtkSMPThreadID tid, T* item )
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

      T* GetLocal( vtkSMPThreadID tid )
        {
        return this->ThreadLocalStorage[tid];
        }

      template<class Derived>
      void FillDerivedThreadLocal ( vtkThreadLocal<Derived>* other )
        {
        for ( typename vtkstd::vector<T*>::size_type i = 0; i < ThreadLocalStorage.size(); ++i )
          {
          other->SetLocal( i, Derived::SafeDownCast(ThreadLocalStorage[i]) );
          }
        }

    protected:
      // the __thread c++Ox modifier cannot be used because this is not static
      // create an explicit map and use the thread id key instead.
      vtkstd::vector<T*> ThreadLocalStorage;
    };

  // ForEach template : parallel loop over an iterator
  void VTK_SMP_EXPORT ForEach(vtkIdType first, vtkIdType last, const vtkFunctor* op );

  void VTK_SMP_EXPORT ForEach(vtkIdType first, vtkIdType last, const vtkFunctorInitialisable* f );

  void VTK_SMP_EXPORT Parallel( const vtkTask* function, const vtkObject* data = NULL, vtkSMPThreadID skipThreads = 1 );

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

}

#endif //__vtkSMP_h__
