/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkMemoryLimitImageDataStreamer.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkMemoryLimitImageDataStreamer.h"

#include "vtkAlgorithmOutput.h"
#include "vtkCommand.h"
#include "vtkExtentTranslator.h"
#include "vtkImageData.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPipelineSize.h"
#include "vtkStreamingDemandDrivenPipeline.h"

vtkCxxRevisionMacro(vtkMemoryLimitImageDataStreamer, "1.9");
vtkStandardNewMacro(vtkMemoryLimitImageDataStreamer);

//----------------------------------------------------------------------------
vtkMemoryLimitImageDataStreamer::vtkMemoryLimitImageDataStreamer()
{
  // Set a default memory limit of 50 Megabytes
  this->MemoryLimit = 50000; 
}


//----------------------------------------------------------------------------
void vtkMemoryLimitImageDataStreamer::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "MemoryLimit (in kb): " << this->MemoryLimit << endl;
}


//----------------------------------------------------------------------------
#ifndef VTK_USE_EXECUTIVES
void vtkMemoryLimitImageDataStreamer::UpdateData(vtkDataObject *out)
{
  // find the right number of pieces
  if (!this->GetInput())
    {
    return;
    }
  
  vtkImageData *input = this->GetInput();
  vtkExtentTranslator *translator = this->GetExtentTranslator();
  translator->SetWholeExtent(out->GetUpdateExtent());

  vtkPipelineSize *sizer = vtkPipelineSize::New();
  this->NumberOfStreamDivisions = 1;
  unsigned long oldSize, size = 0;
  float ratio;
  translator->SetPiece(0);

  // watch for the limiting case where the size is the maximum size
  // represented by an unsigned long. In that case we do not want to do the
  // ratio test. We actual test for size < 0.5 of the max unsigned long which
  // would indicate that oldSize is about at max unsigned long.
  unsigned long maxSize;
  maxSize = (((unsigned long)0x1) << (8*sizeof(unsigned long) - 1));
  
  // we also have to watch how many pieces we are creating. Since
  // NumberOfStreamDivisions is an int, it cannot be more that say 2^31
  // (which is a bit much anyhow) so we also stop if the number of pieces is
  // too large.
  int count = 0;
  
  // double the number of pieces until the size fits in memory
  // or the reduction in size falls to 20%
  do 
    {
    oldSize = size;
    translator->SetNumberOfPieces(this->NumberOfStreamDivisions);
    translator->PieceToExtentByPoints();
    input->SetUpdateExtent(translator->GetExtent());
    input->PropagateUpdateExtent();
    size = sizer->GetEstimatedSize(this->GetInput());
    // watch for the first time through
    if (!oldSize)
      {
      ratio = 0.5;
      }
    // otherwise the normal ratio calculation
    else
      {
      ratio = size/(float)oldSize;
      }
    this->NumberOfStreamDivisions = this->NumberOfStreamDivisions*2;
    count++;
    }
  while (size > this->MemoryLimit && 
         (size < maxSize && ratio < 0.8) && count < 29);
  
  // undo the last *2
  this->NumberOfStreamDivisions = this->NumberOfStreamDivisions/2;
  
  // now call the superclass
  this->vtkImageDataStreamer::UpdateData(out);

  sizer->Delete();
}
#endif

int vtkMemoryLimitImageDataStreamer::ProcessUpstreamRequest(
  vtkInformation *request, 
  vtkInformationVector *inputVector, 
  vtkInformationVector *outputVector)
{
#ifdef VTK_USE_EXECUTIVES
  if(request->Has(vtkStreamingDemandDrivenPipeline::REQUEST_UPDATE_EXTENT()))
    {
    if (this->CurrentDivision == 0)
      {
      // we must set the extent on the input
      vtkInformation* outInfo = outputVector->GetInformationObject(0);
      
      // get the requested update extent
      int outExt[6];
      outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), outExt);
      
      vtkInformation* inInfo = inputVector->GetInformationObject(0)
        ->Get(vtkAlgorithm::INPUT_CONNECTION_INFORMATION())
        ->GetInformationObject(0);
      vtkImageData *input = 
        vtkImageData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
      
      vtkExtentTranslator *translator = this->GetExtentTranslator();
      translator->SetWholeExtent(outExt);
      
      vtkPipelineSize *sizer = vtkPipelineSize::New();
      this->NumberOfStreamDivisions = 1;
      unsigned long oldSize, size = 0;
      float ratio;
      translator->SetPiece(0);
      
      // watch for the limiting case where the size is the maximum size
      // represented by an unsigned long. In that case we do not want to do
      // the ratio test. We actual test for size < 0.5 of the max unsigned
      // long which would indicate that oldSize is about at max unsigned
      // long.
      unsigned long maxSize;
      maxSize = (((unsigned long)0x1) << (8*sizeof(unsigned long) - 1));
      
      // we also have to watch how many pieces we are creating. Since
      // NumberOfStreamDivisions is an int, it cannot be more that say 2^31
      // (which is a bit much anyhow) so we also stop if the number of pieces
      // is too large.
      int count = 0;

      
      // double the number of pieces until the size fits in memory
      // or the reduction in size falls to 20%
      do 
        {
        oldSize = size;
        translator->SetNumberOfPieces(this->NumberOfStreamDivisions);
        translator->PieceToExtentByPoints();

        int inExt[6];
        translator->GetExtent(inExt);
        // set the update extent
        inInfo->Set(
          vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), inExt, 6);
        
        // then propagate it
        vtkAlgorithm *alg = input->GetProducerPort()->GetProducer();
        int index = input->GetProducerPort()->GetIndex();
        vtkStreamingDemandDrivenPipeline *exec = 
          vtkStreamingDemandDrivenPipeline::SafeDownCast(alg->GetExecutive());
        exec->PropagateUpdateExtent(index);

        size = sizer->GetEstimatedSize(this->GetInput());
        // watch for the first time through
        if (!oldSize)
          {
          ratio = 0.5;
          }
        // otherwise the normal ratio calculation
        else
          {
          ratio = size/(float)oldSize;
          }
        this->NumberOfStreamDivisions = this->NumberOfStreamDivisions*2;
        count++;
        }
      while (size > this->MemoryLimit && 
             (size < maxSize && ratio < 0.8) && count < 29);
      
      // undo the last *2
      this->NumberOfStreamDivisions = this->NumberOfStreamDivisions/2;
      }    
    return
      this->Superclass::ProcessUpstreamRequest(request,
                                               inputVector,
                                               outputVector);
    }
  return 0;
#else
  return this->Superclass::ProcessUpstreamRequest(request,
                                                  inputVector,
                                                  outputVector);
#endif
}
