#ifndef __vtkSMPContourFilter_h
#define __vtkSMPContourFilter_h

#include "vtkParallelSMPModule.h" // For export macro
#include "vtkContourFilter.h"
#include "vtkSMPAlgorithm.h"

class VTKPARALLELSMP_EXPORT vtkSMPContourFilter : public vtkContourFilter, public vtkSMPAlgorithm
{
public:
  vtkTypeMacro(vtkSMPContourFilter,vtkContourFilter);
  virtual void PrintSelf(ostream& os, vtkIndent indent);

  static vtkSMPContourFilter *New();

  virtual int ProcessRequest(vtkInformation* request,
                             vtkInformationVector** inInfo,
                             vtkInformationVector* outInfo);
protected:
  vtkSMPContourFilter();
  ~vtkSMPContourFilter();

  virtual int RequestData(vtkInformation* request,
                          vtkInformationVector** inputVector,
                          vtkInformationVector* outputVector);

  virtual int SplitData(vtkInformation* request,
                        vtkInformationVector** inputVector,
                        vtkInformationVector* outputVector,
                        vtkThreadLocal<vtkDataObject>** outputs);

private:
  vtkSMPContourFilter(const vtkSMPContourFilter&);  // Not implemented.
  void operator=(const vtkSMPContourFilter&);  // Not implemented.
};

#endif
