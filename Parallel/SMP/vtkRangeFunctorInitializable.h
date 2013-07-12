#ifndef __vtkRangeFunctorInitializable_h__
#define __vtkRangeFunctorInitializable_h__

#include "vtkParallelSMPModule.h" // For export macro
#include "vtkRangeFunctor.h"
#include <vector> //TODO PIMPLE

class VTKPARALLELSMP_EXPORT vtkRangeFunctorInitializable : public vtkRangeFunctor
{
public:
  vtkTypeMacro(vtkRangeFunctorInitializable,vtkRangeFunctor);
  void PrintSelf(ostream& os, vtkIndent indent);

  virtual void Init ( int ) const = 0;
  bool ShouldInitialize ( int ) const;

protected:
  vtkRangeFunctorInitializable();
  ~vtkRangeFunctorInitializable();

  void Initialized ( int ) const;
  mutable vtkstd::vector<vtkIdType> IsInitialized;

private:
  vtkRangeFunctorInitializable(const vtkRangeFunctorInitializable&);  // Not implemented.
  void operator=(const vtkRangeFunctorInitializable&);  // Not implemented.
};

#endif
