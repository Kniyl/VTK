#ifndef __vtkParallelOperators_h__
#define __vtkParallelOperators_h__

#include "vtkParallelSMPModule.h" // For export macro
#include "vtkObject.h"
#include <vector>

class vtkRangeFunctor;
class vtkRangeFunctorInitializable;
class vtkTreeFunctor;
class vtkTreeFunctorInitializable;
class vtkParallelTree;

class VTKPARALLELSMP_EXPORT vtkParallelOperators : public vtkObject
{
    vtkParallelOperators(const vtkParallelOperators&);
    void operator=(const vtkParallelOperators&);

  protected:
    vtkParallelOperators();
    ~vtkParallelOperators();

  public:
    vtkTypeMacro(vtkParallelOperators,vtkObject);
    static vtkParallelOperators* New();
    void PrintSelf(ostream& os, vtkIndent indent);

    // ForEach template : parallel loop over an iterator
    static void ForEach(
      vtkIdType first, vtkIdType last,
      const vtkRangeFunctor* op, int grain = 0
    );

    static void ForEach(
      vtkIdType first, vtkIdType last,
      const vtkRangeFunctorInitializable* f, int grain = 0
    );
 
    // Same as ForEach but with a guaranteed static partitioning
    // For now it is only an alias, TODO implement it in each
    // runtime.
    static void StaticForEach(
        vtkIdType first, vtkIdType last,
        const vtkRangeFunctor* op, int grain = 0)
      {
      ForEach(first,last,op,grain);
      } 

    static void StaticForEach(
        vtkIdType first, vtkIdType last,
        const vtkRangeFunctorInitializable* op, int grain = 0)
      {
      ForEach(first,last,op,grain);
      } 

    static void Traverse(const vtkParallelTree* Tree, vtkTreeFunctor* func);
};

#endif //__vtkParallelOperators_h__
