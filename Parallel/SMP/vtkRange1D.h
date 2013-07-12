#ifndef __vtkRange1D_h__
#define __vtkRange1D_h__

#include "vtkRange.h"

class VTKPARALLELSMP_EXPORT vtkRange1D : public vtkRange
{
    vtkRange1D(const vtkRange1D&);
    void operator=(const vtkRange1D&);
  protected:
    vtkRange1D();
    ~vtkRange1D();
    
    vtkIdType begin, end;
  public:
    vtkTypeMacro(vtkRange1D,vtkRange);
    static vtkRange1D* New();
    void PrintSelf(ostream& os, vtkIndent indent);

    virtual int GetDimension() const { return 1; }
    void Setup(vtkIdType b, vtkIdType e, int t)
      {
      begin = b;
      end = e;
      tid = t;
      }

    vtkIdType& Begin() { return begin; }
    vtkIdType& End() { return end; }
};

#endif
