#ifndef __vtkTreeFunctorInitializable_h__
#define __vtkTreeFunctorInitializable_h__

#include "vtkParallelSMPModule.h" // For export macro
#include "vtkTreeFunctor.h"
#include <vector> //TODO PIMPLE

class VTKPARALLELSMP_EXPORT vtkTreeFunctorInitializable : public vtkTreeFunctor
{
public:
  vtkTypeMacro(vtkTreeFunctorInitializable,vtkTreeFunctor);
  virtual void PrintSelf(ostream& os, vtkIndent indent);

  virtual void Init ( int ) const = 0;
  bool ShouldInitialize ( int ) const;

protected:
  vtkTreeFunctorInitializable();
  ~vtkTreeFunctorInitializable();

  void Initialized ( int ) const;
  mutable vtkstd::vector<vtkIdType> IsInitialized;

private:
  vtkTreeFunctorInitializable(const vtkTreeFunctorInitializable&);  // Not implemented.
  void operator=(const vtkTreeFunctorInitializable&);  // Not implemented.
};

#endif
