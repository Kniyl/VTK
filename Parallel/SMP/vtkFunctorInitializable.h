#ifndef __vtkFunctorInitializable_h__
#define __vtkFunctorInitializable_h__

#include "vtkParallelSMPModule.h" // For export macro
#include "vtkFunctor.h"
#include <vector> //TODO PIMPLE

class VTKPARALLELSMP_EXPORT vtkFunctorInitializable : public vtkFunctor
{
public:
  vtkTypeMacro(vtkFunctorInitializable,vtkFunctor);
  void PrintSelf(ostream& os, vtkIndent indent);

  virtual void Init ( int ) const = 0;
  bool ShouldInitialize ( int ) const;

protected:
  vtkFunctorInitializable();
  ~vtkFunctorInitializable();

  void Initialized ( int ) const;
  mutable vtkstd::vector<vtkIdType> IsInitialized;

private:
  vtkFunctorInitializable(const vtkFunctorInitializable&);  // Not implemented.
  void operator=(const vtkFunctorInitializable&);  // Not implemented.
};

#endif
