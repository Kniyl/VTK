/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkImageRGBToHSV.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkImageRGBToHSV.h"

#include "vtkImageData.h"
#include "vtkImageProgressIterator.h"
#include "vtkMath.h"
#include "vtkObjectFactory.h"

vtkCxxRevisionMacro(vtkImageRGBToHSV, "1.27");
vtkStandardNewMacro(vtkImageRGBToHSV);

//----------------------------------------------------------------------------
vtkImageRGBToHSV::vtkImageRGBToHSV()
{
  this->Maximum = 255.0;
}

//----------------------------------------------------------------------------
// This templated function executes the filter for any type of data.
template <class T>
void vtkImageRGBToHSVExecute(vtkImageRGBToHSV *self,
                             vtkImageData *inData,
                             vtkImageData *outData,
                             int outExt[6], int id, T *)
{
  vtkImageIterator<T> inIt(inData, outExt);
  vtkImageProgressIterator<T> outIt(outData, outExt, self, id);
  int idxC, maxC;
  float R, G, B, H, S, V;
  float max = self->GetMaximum();
  
  // find the region to loop over
  maxC = inData->GetNumberOfScalarComponents()-1;
  
  // Loop through ouput pixels
  while (!outIt.IsAtEnd())
    {
    T* inSI = inIt.BeginSpan();
    T* outSI = outIt.BeginSpan();
    T* outSIEnd = outIt.EndSpan();
    while (outSI != outSIEnd)
      {
      // Pixel operation
      R = (float)(*inSI) / max; inSI++;
      G = (float)(*inSI) / max; inSI++;
      B = (float)(*inSI) / max; inSI++;

      vtkMath::RGBToHSV(R, G, B, &H, &S, &V);

      H *= max;
      S *= max;
      V *= max;

      if (H > max)
        {
        H = max;
        }
      if (S > max)
        {
        S = max;
        }
      if (V > max)
        {
        V = max;
        }
      
      // assign output.
      *outSI = (T)(H); outSI++;
      *outSI = (T)(S); outSI++;
      *outSI = (T)(V); outSI++;
      
      for (idxC = 3; idxC <= maxC; idxC++)
        {
        *outSI++ = *inSI++;
        }
      }
    inIt.NextSpan();
    outIt.NextSpan();
    }
}

//----------------------------------------------------------------------------
void vtkImageRGBToHSV::ThreadedExecute(vtkImageData *inData, 
                                         vtkImageData *outData,
                                         int outExt[6], int id)
{
  vtkDebugMacro(<< "Execute: inData = " << inData 
  << ", outData = " << outData);
  
  // this filter expects that input is the same type as output.
  if (inData->GetScalarType() != outData->GetScalarType())
    {
    vtkErrorMacro(<< "Execute: input ScalarType, " << inData->GetScalarType()
    << ", must match out ScalarType " << outData->GetScalarType());
    return;
    }
  
  // need three components for input and output
  if (inData->GetNumberOfScalarComponents() < 3)
    {
    vtkErrorMacro("Input has too few components");
    return;
    }
  if (outData->GetNumberOfScalarComponents() < 3)
    {
    vtkErrorMacro("Output has too few components");
    return;
    }

  switch (inData->GetScalarType())
    {
    vtkTemplateMacro6(vtkImageRGBToHSVExecute, this, inData, 
                      outData, outExt, id, static_cast<VTK_TT *>(0));
    default:
      vtkErrorMacro(<< "Execute: Unknown ScalarType");
      return;
    }
}

void vtkImageRGBToHSV::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "Maximum: " << this->Maximum << "\n";
}

