/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkImageLuminance.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkImageLuminance.h"

#include "vtkImageData.h"
#include "vtkImageProgressIterator.h"
#include "vtkObjectFactory.h"

#include <math.h>

vtkCxxRevisionMacro(vtkImageLuminance, "1.22.10.1");
vtkStandardNewMacro(vtkImageLuminance);

//----------------------------------------------------------------------------
// This method overrides information set by parent's ExecuteInformation.
void vtkImageLuminance::ExecuteInformation(vtkImageData *vtkNotUsed(inData), 
                                           vtkImageData *outData)
{
  outData[0]->SetNumberOfScalarComponents(1);
}

//----------------------------------------------------------------------------
// This execute method handles boundaries.
// it handles boundaries. Pixels are just replicated to get values 
// out of extent.
template <class T>
void vtkImageLuminanceExecute(vtkImageLuminance *self, vtkImageData *inData,
                              vtkImageData *outData,
                              int outExt[6], int id, T *)
{
  vtkImageIterator<T> inIt(inData, outExt);
  vtkImageProgressIterator<T> outIt(outData, outExt, self, id);
  float luminance;

  // Loop through ouput pixels
  while (!outIt.IsAtEnd())
    {
    T* inSI = inIt.BeginSpan();
    T* outSI = outIt.BeginSpan();
    T* outSIEnd = outIt.EndSpan();
    while (outSI != outSIEnd)
      {
      // now process the components
      luminance =  0.30 * *inSI++;
      luminance += 0.59 * *inSI++;
      luminance += 0.11 * *inSI++;
      *outSI = static_cast<T>(luminance);
      ++outSI;
      }
    inIt.NextSpan();
    outIt.NextSpan();
    }
}


//----------------------------------------------------------------------------
// This method contains a switch statement that calls the correct
// templated function for the input data type.  The output data
// must match input type.  This method does handle boundary conditions.
void vtkImageLuminance::ThreadedExecute (vtkImageData ***inData, 
                                        vtkImageData **outData,
                                        int outExt[6], int id)
{
  vtkDebugMacro(<< "Execute: inData = " << inData 
  << ", outData = " << outData);
  
  // this filter expects that input is the same type as output.
  if (inData[0][0]->GetNumberOfScalarComponents() != 3)
    {
    vtkErrorMacro(<< "Execute: input must have 3 components, but has " << inData[0][0]->GetNumberOfScalarComponents());
    return;
    }
  
  // this filter expects that input is the same type as output.
  if (inData[0][0]->GetScalarType() != outData[0]->GetScalarType())
    {
    vtkErrorMacro(<< "Execute: input ScalarType, " << inData[0][0]->GetScalarType()
    << ", must match out ScalarType " << outData[0]->GetScalarType());
    return;
    }
  
  switch (inData[0][0]->GetScalarType())
    {
    vtkTemplateMacro6(vtkImageLuminanceExecute, this, inData[0][0], outData[0], 
                      outExt, id, static_cast<VTK_TT *>(0));
    default:
      vtkErrorMacro(<< "Execute: Unknown ScalarType");
      return;
    }
}










