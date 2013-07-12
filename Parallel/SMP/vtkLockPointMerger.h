#ifndef _vtkLockPointMerger_h_
#define _vtkLockPointMerger_h_

//extra synchonization locks for case when you lack allocators per CPU

#include "vtkParallelSMPModule.h" // For export macro
#include "vtkRangeFunctor.h"
#include "vtkRange1D.h"

class vtkPoints;
class vtkPointData;
class vtkDummyMergeFunctor;

struct VTKPARALLELSMP_EXPORT vtkLockPointMerger : public vtkRangeFunctor
{
  vtkDummyMergeFunctor* Functor;
  vtkIdType NumberOfPointsFirstThread;

  vtkTypeMacro(vtkLockPointMerger,vtkRangeFunctor);
  static vtkLockPointMerger* New();
  void PrintSelf(ostream &os, vtkIndent indent);

  void operator()(vtkRange*) const;

protected:
  vtkLockPointMerger() {}
  ~vtkLockPointMerger() {}

private:
  vtkLockPointMerger(const vtkLockPointMerger&);
  void operator =(const vtkLockPointMerger&);
};

#endif //_vtkLockPointMerger_h_
