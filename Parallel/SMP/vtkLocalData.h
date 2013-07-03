#ifndef __vtkLocalData_h__
#define __vtkLocalData_h__

#include "vtkParallelSMPModule.h" // For export macro
#include "vtkObject.h"

class VTKPARALLELSMP_EXPORT vtkLocalData : public vtkObject
{
    vtkLocalData(const vtkLocalData&);
    void operator=(const vtkLocalData&);

  protected:
    vtkLocalData();
    ~vtkLocalData();

  public:
    vtkTypeMacro(vtkLocalData,vtkObject);
    void PrintSelf(ostream& os, vtkIndent indent);
    static vtkLocalData* New();
};

#endif
