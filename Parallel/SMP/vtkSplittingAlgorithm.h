#ifndef __vtkSplittingAlgorithm_h
#define __vtkSplittingAlgorithm_h

#include "vtkParallelSMPModule.h" // For export macro
#include "vtkSetGet.h"

class vtkDataObject;
class vtkInformation;
class vtkInformationVector;
template<class T> class vtkThreadLocal;

#define vtkSplittingFilterStandardProcessRequest(classname)               \
int classname::ProcessRequest(vtkInformation* request,                    \
                                        vtkInformationVector** inVector,  \
                                        vtkInformationVector* outVector)  \
  {                                                                       \
  if(!this->vtkSplittingAlgorithm::ProcessRequest(                        \
        request,inVector,outVector,this->GetNumberOfOutputPorts()))       \
    {                                                                     \
    return this->Superclass::ProcessRequest(request,inVector,outVector);  \
    }                                                                     \
  return 1;                                                               \
  }


class VTKPARALLELSMP_EXPORT vtkSplittingAlgorithm
{
public:
  virtual int ProcessRequest(vtkInformation* request,
                             vtkInformationVector** inInfo,
                             vtkInformationVector* outInfo,
                             int NumberOfOutputPorts);

  vtkBooleanMacro(SplitDataset,int);
  void SetSplitDataset(int value) { this->SplitDataset = value; }
  int GetSplitDataset() { return this->SplitDataset; }

protected:
  vtkSplittingAlgorithm();
  ~vtkSplittingAlgorithm();

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
  vtkSplittingAlgorithm(const vtkSplittingAlgorithm&);  // Not implemented.
  void operator=(const vtkSplittingAlgorithm&);  // Not implemented.
};

#endif
