#ifndef _vtkParallelPointMerger_h_
#define _vtkParallelPointMerger_h_

#include "vtkParallelSMPModule.h" // For export macro
#include "vtkTask.h"

class vtkDummyMergeFunctor;
class vtkSMPMergePoints;

struct VTKPARALLELSMP_EXPORT vtkParallelPointMerger : public vtkTask
{

  vtkTypeMacro(vtkParallelPointMerger,vtkTask);
  static vtkParallelPointMerger* New();
  void PrintSelf(ostream &os, vtkIndent indent);
  void Execute( vtkSMPMergePoints* locator ) const;
  
  void SetUsefullData(vtkDummyMergeFunctor* f, vtkIdType** t);

protected:
  vtkParallelPointMerger();
  ~vtkParallelPointMerger();

  int MustTreatBucket( vtkIdType idx ) const;
  vtkIdType** TreatedTable;
  vtkDummyMergeFunctor* self;

private:
  vtkParallelPointMerger(const vtkParallelPointMerger&);
  void operator =(const vtkParallelPointMerger&);
};

#endif //_vtkParallelPointMerger_h_
