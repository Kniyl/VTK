#ifndef __vtkSMPClipDataSet_h
#define __vtkSMPClipDataSet_h

#include "vtkParallelSMPModule.h" // For export macro
#include "vtkClipDataSet.h"
#include "vtkSplittingAlgorithm.h"

class vtkDataObject;
template<class T> class vtkThreadLocal;

class VTKPARALLELSMP_EXPORT vtkSMPClipDataSet : public vtkClipDataSet, public vtkSplittingAlgorithm
{
public:
  vtkTypeMacro(vtkSMPClipDataSet,vtkClipDataSet);
  void PrintSelf(ostream& os, vtkIndent indent);

  static vtkSMPClipDataSet *New();
  virtual int ProcessRequest(vtkInformation* request,
                             vtkInformationVector** inVec,
                             vtkInformationVector* outVec);

protected:
  vtkSMPClipDataSet();
  ~vtkSMPClipDataSet();

  virtual int RequestData(vtkInformation* request,
                          vtkInformationVector** inputVector,
                          vtkInformationVector* outputVector);

  virtual int SplitData(vtkInformation* request,
                        vtkInformationVector** inputVector,
                        vtkInformationVector* outputVector,
                        vtkThreadLocal<vtkDataObject>** outputData);

private:
  vtkSMPClipDataSet(const vtkSMPClipDataSet&);  // Not implemented.
  void operator=(const vtkSMPClipDataSet&);  // Not implemented.
};

#endif
