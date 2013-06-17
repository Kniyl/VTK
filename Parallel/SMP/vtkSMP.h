#ifndef __vtkSMP_h__
#define __vtkSMP_h__

#include "vtkParallelSMPModule.h" // For export macro
#include "vtkObject.h"
#include <vector>

#include "vtkFunctorInitializable.h"

class vtkCellArray;
class vtkCellData;
class vtkFunctor;
class vtkFunctorIntializable;
class vtkParallelTree;
class vtkPointData;
class vtkPoints;
class vtkSMPMergePoints;
class vtkTask;

//======================================================================================
namespace vtkSMP
{
  int InternalGetNumberOfThreads();
  int InternalGetTid();

  //======================================================================================
  template<class T>
  class VTKPARALLELSMP_EXPORT vtkThreadLocal : public vtkObject
  {
    protected :
      vtkThreadLocal() : vtkObject(), ThreadLocalStorage(InternalGetNumberOfThreads(), NULL) {}
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
      typedef typename vtkstd::vector<T*>::iterator iterator;

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

      T* NewLocal ( T* specificImpl )
        {
        int tid = InternalGetTid();
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

      T* NewLocal ( )
        {
        int tid = InternalGetTid();
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
      void FillDerivedThreadLocal( vtkThreadLocal<Derived>* other )
        {
        T* elem;
        iterator src = ThreadLocalStorage.begin();
        for ( typename vtkSMP::vtkThreadLocal<Derived>::iterator it = other->Begin();
              it != other->End(); ++it, ++src )
          {
          if ( (elem = *it) ) elem->UnRegister(other);
          Derived* d = (*it) = Derived::SafeDownCast(*src);
          if ( d ) d->Register(other);
          }
        }

    protected:
      vtkstd::vector<T*> ThreadLocalStorage;
    };

  //======================================================================================
  // ForEach template : parallel loop over an iterator
  void VTKPARALLELSMP_EXPORT ForEach( vtkIdType first, vtkIdType last, const vtkFunctor* op, int grain = 0 );

  void VTKPARALLELSMP_EXPORT ForEach( vtkIdType first, vtkIdType last, const vtkFunctorInitializable* f, int grain = 0 );

  template<class T>
  void VTKPARALLELSMP_EXPORT Parallel( const vtkTask* function,
                                typename vtkSMP::vtkThreadLocal<T>::iterator data1,
                                vtkIdType skipThreads = 1 );

  template<class T1, class T2, class T3, class T4, class T5, class T6>
  void VTKPARALLELSMP_EXPORT Parallel( const vtkTask* function,
                                typename vtkSMP::vtkThreadLocal<T1>::iterator data1,
                                typename vtkSMP::vtkThreadLocal<T2>::iterator data2,
                                typename vtkSMP::vtkThreadLocal<T3>::iterator data3,
                                typename vtkSMP::vtkThreadLocal<T4>::iterator data4,
                                typename vtkSMP::vtkThreadLocal<T5>::iterator data5,
                                typename vtkSMP::vtkThreadLocal<T6>::iterator data6,
                                vtkstd::vector<vtkIdType>::iterator offset1,
                                vtkstd::vector<vtkIdType>::iterator offset2,
                                vtkstd::vector<vtkIdType>::iterator offset3,
                                vtkstd::vector<vtkIdType>::iterator offset4,
                                vtkstd::vector<vtkIdType>::iterator offset5,
                                vtkstd::vector<vtkIdType>::iterator offset6,
                                vtkstd::vector<vtkIdType>::iterator offset7,
                                vtkstd::vector<vtkIdType>::iterator offset8,
                                vtkIdType skipThreads = 1 );

  void VTKPARALLELSMP_EXPORT Traverse( const vtkParallelTree* Tree, vtkFunctor* func );

  void VTKPARALLELSMP_EXPORT MergePoints( vtkPoints* outPoints, vtkThreadLocal<vtkPoints>* inPoints, const double bounds[6],
                                   vtkPointData* outPtsData, vtkThreadLocal<vtkPointData>* inPtsData,
                                   vtkCellArray* outVerts, vtkThreadLocal<vtkCellArray>* inVerts,
                                   vtkCellArray* outLines, vtkThreadLocal<vtkCellArray>* inLines,
                                   vtkCellArray* outPolys, vtkThreadLocal<vtkCellArray>* inPolys,
                                   vtkCellArray* outStrips, vtkThreadLocal<vtkCellArray>* inStrips,
                                   vtkCellData* outCellsData, vtkThreadLocal<vtkCellData>* inCellsData, int SkipThreads );

  void VTKPARALLELSMP_EXPORT MergePoints( vtkSMPMergePoints* outPoints, vtkThreadLocal<vtkSMPMergePoints>* inPoints,
                                   vtkPointData* outPtsData, vtkThreadLocal<vtkPointData>* inPtsData,
                                   vtkCellArray* outVerts, vtkThreadLocal<vtkCellArray>* inVerts,
                                   vtkCellArray* outLines, vtkThreadLocal<vtkCellArray>* inLines,
                                   vtkCellArray* outPolys, vtkThreadLocal<vtkCellArray>* inPolys,
                                   vtkCellArray* outStrips, vtkThreadLocal<vtkCellArray>* inStrips,
                                   vtkCellData* outCellsData, vtkThreadLocal<vtkCellData>* inCellsData, int SkipThreads );

}

#endif //__vtkSMP_h__
