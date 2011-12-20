#ifndef __vtkSMPContourFilter_h
#define __vtkSMPContourFilter_h

#include "vtkContourFilter.h"

class VTK_SMP_EXPORT vtkSMPContourFilter : public vtkContourFilter
{
public:
  vtkTypeMacro(vtkSMPContourFilter,vtkContourFilter);
  void PrintSelf(ostream& os, vtkIndent indent);

  static vtkSMPContourFilter *New();

protected:
  vtkSMPContourFilter();
  ~vtkSMPContourFilter();

  virtual int RequestData(vtkInformation* request,
                          vtkInformationVector** inputVector,
                          vtkInformationVector* outputVector);

private:
  vtkSMPContourFilter(const vtkSMPContourFilter&);  // Not implemented.
  void operator=(const vtkSMPContourFilter&);  // Not implemented.
};

#endif
