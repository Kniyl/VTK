#ifndef __vtkFunctor_h__
#define __vtkFunctor_h__

#include "vtkParallelSMPModule.h" // For export macro
#include "vtkLocalData.h"

class VTKPARALLELSMP_EXPORT vtkFunctor : public vtkObject
{
  vtkFunctor(const vtkFunctor&);  // Not implemented.
  void operator=(const vtkFunctor&);  // Not implemented.

  void ComputeMasterTID();

protected:
  vtkFunctor();
  ~vtkFunctor();

  int MasterThreadId;

public:
  vtkTypeMacro(vtkFunctor,vtkObject);
  void PrintSelf(ostream& os, vtkIndent indent);

  virtual void operator()(vtkIdType, vtkLocalData*) const = 0;
  virtual vtkLocalData* getLocal(int tid) const;
};

#endif
