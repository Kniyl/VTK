#include "vtkSplittingAlgorithm.h"
#include "vtkParallelCompositeDataPipeline.h"
#include "vtkThreadLocal.h"

#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkDataObject.h"
#include "vtkMultiPieceDataSet.h"


//----------------------------------------------------------------------------
vtkSplittingAlgorithm::vtkSplittingAlgorithm()
  {
  SplitDataset = 0;
  }

//----------------------------------------------------------------------------
vtkSplittingAlgorithm::~vtkSplittingAlgorithm()
  {
  }

//----------------------------------------------------------------------------
// This is the superclasses style of Execute method.  Convert it into
// an SMP style Execute method.
int vtkSplittingAlgorithm::ProcessRequest(
  vtkInformation* request,
  vtkInformationVector** inputVector,
  vtkInformationVector* outputVector,
  int NumberOfOutputPorts)
  {
  if( !request->Has(vtkDemandDrivenPipeline::REQUEST_DATA()) ||
      !this->SplitDataset ||
      !outputVector->GetInformationObject(0)->Has(
        vtkParallelCompositeDataPipeline::DATA_OBJECT_CONCRETE_TYPE()) )
    {
    return 0;
    }
  vtkThreadLocal<vtkDataObject>** outputs =
        new vtkThreadLocal<vtkDataObject>*[NumberOfOutputPorts];
  int i;
  for (i = 0; i < NumberOfOutputPorts; ++i)
    {
    outputs[i] = vtkThreadLocal<vtkDataObject>::New();
    vtkInformation* outInfo = outputVector->GetInformationObject(i);
    outputs[i]->SetSpecificClassName(
      outInfo->Get(vtkParallelCompositeDataPipeline::DATA_OBJECT_CONCRETE_TYPE()));
    }

  if (this->SplitData(request, inputVector, outputVector, outputs))
    {
    for (i = 0; i < NumberOfOutputPorts; ++i)
      {
      vtkMultiPieceDataSet* ds = vtkMultiPieceDataSet::SafeDownCast(
            outputVector->GetInformationObject(i)->Get(
              vtkDataObject::DATA_OBJECT()));
      ds->SetNumberOfPieces(1);
      unsigned int pieceNum = 0;
      for (vtkThreadLocal<vtkDataObject>::iterator outIter = outputs[i]->Begin();
           outIter != outputs[i]->End(); ++outIter)
        {
        ds->SetPiece(pieceNum++, *outIter);
        }

      outputs[i]->FastDelete();
      }
    delete [] outputs;

    return 1;
    }

  return 0;
}

//----------------------------------------------------------------------------
// The execute method created by the subclass.
int vtkSplittingAlgorithm::SplitData(
  vtkInformation* vtkNotUsed( request ),
  vtkInformationVector** vtkNotUsed( inputVector ),
  vtkInformationVector* vtkNotUsed( outputVector ),
  vtkThreadLocal<vtkDataObject>** vtkNotUsed( outputData ))
{
  return 1;
}
