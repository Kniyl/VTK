/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkImageShiftScale.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkImageShiftScale.h"

#include "vtkImageData.h"
#include "vtkImageProgressIterator.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"

vtkCxxRevisionMacro(vtkImageShiftScale, "1.48");
vtkStandardNewMacro(vtkImageShiftScale);

//----------------------------------------------------------------------------
// Constructor sets default values
vtkImageShiftScale::vtkImageShiftScale()
{
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
  this->Shift = 0.0;
  this->Scale = 1.0;
  this->OutputScalarType = -1;
  this->ClampOverflow = 0;
}



//----------------------------------------------------------------------------
void vtkImageShiftScale::ExecuteInformation (
  vtkInformation * vtkNotUsed(request),
  vtkInformationVector ** vtkNotUsed( inputVector ), 
  vtkInformationVector * outputVector)
{
  if (this->OutputScalarType != -1)
    {
    // get the info objects
    vtkInformation* outInfo = outputVector->GetInformationObject(0);
    outInfo->Set(vtkDataObject::SCALAR_TYPE(),this->OutputScalarType);
    }
}




//----------------------------------------------------------------------------
// This templated function executes the filter for any type of data.
template <class IT, class OT>
void vtkImageShiftScaleExecute(vtkImageShiftScale *self,
                               vtkImageData *inData,
                               vtkImageData *outData,
                               int outExt[6], int id,  IT *, OT *)
{
  vtkImageIterator<IT> inIt(inData, outExt);
  vtkImageProgressIterator<OT> outIt(outData, outExt, self, id);
  double typeMin, typeMax, val;
  int clamp;
  double shift = self->GetShift();
  double scale = self->GetScale();

  // for preventing overflow
  typeMin = outData->GetScalarTypeMin();
  typeMax = outData->GetScalarTypeMax();
  clamp = self->GetClampOverflow();
    
  // Loop through ouput pixels
  while (!outIt.IsAtEnd())
    {
    IT* inSI = inIt.BeginSpan();
    OT* outSI = outIt.BeginSpan();
    OT* outSIEnd = outIt.EndSpan();
    if (clamp)
      {
      while (outSI != outSIEnd)
        {
        // Pixel operation
        val = ((double)(*inSI) + shift) * scale;
        if (val > typeMax)
          {
          val = typeMax;
          }
        if (val < typeMin)
          {
          val = typeMin;
          }
        *outSI = (OT)(val);
        ++outSI;
        ++inSI;
        }
      }
    else
      {
      while (outSI != outSIEnd)
        {
        // Pixel operation
        *outSI = (OT)(((double)(*inSI) + shift) * scale);
        ++outSI;
        ++inSI;
        }
      }
    inIt.NextSpan();
    outIt.NextSpan();
    }
}



//----------------------------------------------------------------------------
template <class T>
void vtkImageShiftScaleExecute1(vtkImageShiftScale *self,
                                vtkImageData *inData,
                                vtkImageData *outData,
                                int outExt[6], int id, T *)
{
  switch (outData->GetScalarType())
    {
    vtkTemplateMacro7(vtkImageShiftScaleExecute, self, inData,
                      outData,outExt, id,
                      static_cast<T *>(0), static_cast<VTK_TT *>(0));
    default:
      vtkGenericWarningMacro("Execute: Unknown input ScalarType");
      return;
    }
}




//----------------------------------------------------------------------------
// This method is passed a input and output data, and executes the filter
// algorithm to fill the output from the input.
// It just executes a switch statement to call the correct function for
// the datas data types.
void vtkImageShiftScale::ThreadedExecute (vtkImageData *inData, 
                                         vtkImageData *outData,
                                         int outExt[6], int id)
{
  switch (inData->GetScalarType())
    {
    vtkTemplateMacro6(vtkImageShiftScaleExecute1, this, 
                      inData, outData, outExt, id, static_cast<VTK_TT *>(0));
    default:
      vtkErrorMacro(<< "Execute: Unknown ScalarType");
      return;
    }
}



void vtkImageShiftScale::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "Shift: " << this->Shift << "\n";
  os << indent << "Scale: " << this->Scale << "\n";
  os << indent << "Output Scalar Type: " << this->OutputScalarType << "\n";
  os << indent << "ClampOverflow: ";
  if (this->ClampOverflow)
    {
    os << "On\n";
    }
  else 
    {
    os << "Off\n";
    }
}

