#ifndef __vtkSMPAlgorithm_h
#define __vtkSMPAlgorithm_h

#include "vtkParallelSMPModule.h" // For export macro
#include "vtkSetGet.h"

class vtkDataObject;
class vtkInformation;
class vtkInformationVector;
template<class T> class vtkThreadLocal;

class VTKPARALLELSMP_EXPORT vtkSMPAlgorithm
{
public:
  virtual int ProcessRequest(vtkInformation* request,
                             vtkInformationVector** inInfo,
                             vtkInformationVector* outInfo,
                             int NumberOfOutputPorts);

  vtkBooleanMacro(SplitDataset,int);
  void SetSplitDataset(int value) { this->SplitDataset; }
  int GetSplitDataset() { return this->SplitDataset; }

protected:
  vtkSMPAlgorithm();
  ~vtkSMPAlgorithm();

  // Description:
  // This is called by the superclass.
  // This is the method you should override.
  // Use the inputVector to populate the outputData, the
  // superclass will be responsible for merging each local
  // output into the outputVector.
  virtual int SplitData(vtkInformation* request,
                        vtkInformationVector** inputVector,
                        vtkInformationVector* outputVector,
                        vtkThreadLocal<vtkDataObject>** outputData);

  int SplitDataset;

private:
  vtkSMPAlgorithm(const vtkSMPAlgorithm&);  // Not implemented.
  void operator=(const vtkSMPAlgorithm&);  // Not implemented.
};

#endif
