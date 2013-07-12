#ifndef __vtkTreeFunctor_h__
#define __vtkTreeFunctor_h__

#include "vtkParallelSMPModule.h" // For export macro
#include "vtkLocalData.h"

class VTKPARALLELSMP_EXPORT vtkTreeFunctor : public vtkObject
{
  vtkTreeFunctor(const vtkTreeFunctor&);  // Not implemented.
  void operator=(const vtkTreeFunctor&);  // Not implemented.

  void ComputeMasterTID();

protected:
  vtkTreeFunctor();
  ~vtkTreeFunctor();

  int MasterThreadId;

public:
  vtkTypeMacro(vtkTreeFunctor,vtkObject);
  void PrintSelf(ostream& os, vtkIndent indent);

  virtual void operator()(vtkIdType, vtkLocalData*) const = 0;
  virtual vtkLocalData* getLocal(int tid) const;
};

#endif
