#ifndef __vtkRange_h__
#define __vtkRange_h__

#include "vtkObject.h"
#include "vtkParallelSMPModule.h" // For export macro

class VTKPARALLELSMP_EXPORT vtkRange : public vtkObject
{
    vtkRange(const vtkRange&);
    void operator=(const vtkRange&);
  protected:
    vtkRange();
    ~vtkRange();

    int tid;
  public:
    vtkTypeMacro(vtkRange,vtkObject);
    void PrintSelf(ostream& os, vtkIndent indent);

    virtual int GetDimension() const = 0;

    void SetTid(int t) { tid = t; }
    int& GetTid() { return tid; }
};

#endif
