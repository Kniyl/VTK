#include "vtkSMPAlgorithm.h"
#include "vtkSMPPipeline.h"
#include "vtkThreadLocal.h"

#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkDataObject.h"
#include "vtkMultiPieceDataSet.h"


//----------------------------------------------------------------------------
vtkSMPAlgorithm::vtkSMPAlgorithm()
  {
  SplitDataset = 0;
  }

//----------------------------------------------------------------------------
vtkSMPAlgorithm::~vtkSMPAlgorithm()
  {
  }

//----------------------------------------------------------------------------
// This is the superclasses style of Execute method.  Convert it into
// an SMP style Execute method.
int vtkSMPAlgorithm::ProcessRequest(
  vtkInformation* request,
  vtkInformationVector** inputVector,
  vtkInformationVector* outputVector,
  int NumberOfOutputPorts)
  {
  if(!request->Has(vtkDemandDrivenPipeline::REQUEST_DATA()) || !this->SplitDataset)
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
    if (outInfo->Has(vtkSMPPipeline::DATA_OBJECT_CONCRETE_TYPE()))
      {
      outputs[i]->SetSpecificClassName(
          outInfo->Get(vtkSMPPipeline::DATA_OBJECT_CONCRETE_TYPE()));
      }
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
int vtkSMPAlgorithm::SplitData(
  vtkInformation* vtkNotUsed( request ),
  vtkInformationVector** vtkNotUsed( inputVector ),
  vtkInformationVector* vtkNotUsed( outputVector ),
  vtkThreadLocal<vtkDataObject>** vtkNotUsed( outputData ))
{
  return 1;
}
