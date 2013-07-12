#ifndef __vtkRangeFunctor_h__
#define __vtkRangeFunctor_h__

#include "vtkParallelSMPModule.h" // For export macro
#include "vtkRange.h"

class VTKPARALLELSMP_EXPORT vtkRangeFunctor : public vtkObject
{
  vtkRangeFunctor(const vtkRangeFunctor&);  // Not implemented.
  void operator=(const vtkRangeFunctor&);  // Not implemented.

  void ComputeMasterTID();

protected:
  vtkRangeFunctor();
  ~vtkRangeFunctor();

  int MasterThreadId;

public:
  vtkTypeMacro(vtkRangeFunctor,vtkObject);
  void PrintSelf(ostream& os, vtkIndent indent);

  virtual void operator()(vtkRange*) const = 0;
};

#endif
