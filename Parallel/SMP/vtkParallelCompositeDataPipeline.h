#ifndef VTKSMPPIPELINE_H
#define VTKSMPPIPELINE_H

#include "vtkParallelSMPModule.h" // For export macro
#include "vtkCompositeDataPipeline.h"

class VTKPARALLELSMP_EXPORT vtkParallelCompositeDataPipeline : public vtkCompositeDataPipeline
{
  vtkParallelCompositeDataPipeline(const vtkParallelCompositeDataPipeline&);
  void operator =(const vtkParallelCompositeDataPipeline&);

protected:
  vtkParallelCompositeDataPipeline();
  ~vtkParallelCompositeDataPipeline();

  virtual int ExecuteData(vtkInformation* request,
                          vtkInformationVector** inInfoVec,
                          vtkInformationVector* outInfoVec);
  virtual void ExecuteSimpleAlgorithm(vtkInformation* request,
                                      vtkInformationVector** inInfoVec,
                                      vtkInformationVector* outInfoVec,
                                      int compositePort);
  // Check whether the data object in the pipeline information for an
  // output port exists and has a valid type.
  virtual int CheckDataObject(int port, vtkInformationVector* outInfo);

public:
  friend class ParallelFilterExecutor;

  vtkTypeMacro(vtkParallelCompositeDataPipeline, vtkCompositeDataPipeline);
  static vtkParallelCompositeDataPipeline* New();
  void PrintSelf(ostream &os, vtkIndent indent);

  // Description:
  // Key defining the concrete type of output data to make sure
  // the vtkSplittingAlgorithm will produce the right type of temp data.
  static vtkInformationStringKey* DATA_OBJECT_CONCRETE_TYPE();

};

#endif // VTKSMPPIPELINE_H
