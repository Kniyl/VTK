/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkImageRectilinearWipe.cxx
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

  Copyright (c) 1993-2002 Ken Martin, Will Schroeder, Bill Lorensen 
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkImageRectilinearWipe.h"

#include "vtkImageData.h"
#include "vtkObjectFactory.h"

vtkCxxRevisionMacro(vtkImageRectilinearWipe, "1.1.4.1");
vtkStandardNewMacro(vtkImageRectilinearWipe);

//----------------------------------------------------------------------------
vtkImageRectilinearWipe::vtkImageRectilinearWipe()
{
  this->Position[0] = 0;
  this->Position[1] = 0;
  this->Wipe = VTK_WIPE_QUAD;
}

//----------------------------------------------------------------------------
// This templated function executes the filter for any type of data.
// Handles the two input operations
template <class T>
void vtkImageRectilinearWipeExecute2(vtkImageRectilinearWipe *self,
                           vtkImageData *inData, T *inPtr,
                           vtkImageData *outData, 
                           T *outPtr,
                           int outExt[6], int id)
{
  int idxR, idxY, idxZ;
  int maxY, maxZ;
  int inIncX, inIncY, inIncZ;
  int outIncX, outIncY, outIncZ;
  int rowLength;
  unsigned long count = 0;
  unsigned long target;

  // find the region to loop over
  rowLength = (outExt[1] - outExt[0]+1)*inData->GetNumberOfScalarComponents();
  maxY = outExt[3] - outExt[2]; 
  maxZ = outExt[5] - outExt[4];
    
  target = (unsigned long)((maxZ+1)*(maxY+1)/50.0);
  target++;
  
  // Get increments to march through data 
  inData->GetContinuousIncrements(outExt, inIncX, inIncY, inIncZ);
  outData->GetContinuousIncrements(outExt, outIncX, outIncY, outIncZ);

  // Loop through output pixels
  for (idxZ = 0; idxZ <= maxZ; idxZ++)
    {
    for (idxY = 0; idxY <= maxY; idxY++)
      {
      if (!id) 
        {
        if (!(count%target))
          {
          self->UpdateProgress(count/(50.0*target));
          }
        count++;
        }
      for (idxR = 0; idxR < rowLength; idxR++)
        {
        *outPtr = *inPtr;
        outPtr++;
        inPtr++;
        }
      outPtr += outIncY;
      inPtr += inIncY;
      }
    outPtr += outIncZ;
    inPtr += inIncZ;
    }
}



//----------------------------------------------------------------------------
// This function adjusts the extents of the wipe to the output extents.
int vtkImageRectilinearWipeClampExtents(int wipeExt[6], int outExt[6])
{
  int status = 1;
  
  for (int i = 0; i < 3; i++)
    {
    // the lower and upper extents cannot be below the lower output extent
    if (wipeExt[2*i] < outExt[2*i])
      {
      wipeExt[2*i] = outExt[2*i];
      }
    if (wipeExt[2*i + 1] < outExt[2*i])
      {
      wipeExt[2*i + 1] = outExt[2*i];
      status = 0;
      }

    // the lower and upper extents cannot be above the upper output extent
    if (wipeExt[2*i] > outExt[2*i + 1])
      {
      wipeExt[2*i] = outExt[2*i + 1];
      status = 0;
      }
    if (wipeExt[2*i + 1] > outExt[2*i + 1])
      {
      wipeExt[2*i + 1] = outExt[2*i + 1];
      }
    }
  return status;
}
//----------------------------------------------------------------------------
// This method is passed a input and output regions, and executes the filter
// algorithm to fill the output from the inputs based on the Wipe ivar.
void vtkImageRectilinearWipe::ThreadedExecute(vtkImageData **inData, 
                                           vtkImageData *outData,
                                           int outExt[6], int id)
{
  void *inPtr;
  void *outPtr;
  int wipeExt[6];
  int wholeExt[6];
  int whichInput = 0;
  
  vtkDebugMacro(<< "Execute: inData = " << inData 
                << ", outData = " << outData);
  
  // Make sure the inputs/output are valid
  if (inData[0] == NULL)
    {
    vtkErrorMacro(<< "Input " << 0 << " must be specified.");
    return;
    }
  
  // this filter expects that input is the same type as output.
  if (inData[0]->GetScalarType() != outData->GetScalarType())
    {
    vtkErrorMacro(<< "Execute: input ScalarType, "
                  << inData[0]->GetScalarType()
                  << ", must match out ScalarType "
                  << outData->GetScalarType());
    return;
    }

  if (inData[1] == NULL)
    {
    vtkErrorMacro(<< "Input " << 1 << " must be specified.");
    return;
    }

  // this filter expects that inputs that have the same number of components
  if (inData[0]->GetNumberOfScalarComponents() != 
      inData[1]->GetNumberOfScalarComponents())
    {
    vtkErrorMacro(<< "Execute: input1 NumberOfScalarComponents, "
                  << inData[0]->GetNumberOfScalarComponents()
                  << ", must match out input2 NumberOfScalarComponents "
                  << inData[1]->GetNumberOfScalarComponents());
    return;
    }
  
  // Wipe pattern depends on the whole extent.
  outData->GetWholeExtent(wholeExt);

  // Each quadrant is processed separately
  // lower left
  wipeExt[0] = wholeExt[0];
  wipeExt[1] = wholeExt[0] + this->Position[0];
  wipeExt[2] = wholeExt[2];
  wipeExt[3] = wholeExt[2] + this->Position[1];
  wipeExt[4] = wholeExt[4];
  wipeExt[5] = wholeExt[5];
  if (vtkImageRectilinearWipeClampExtents(wipeExt, outExt))
    {

    outPtr = outData->GetScalarPointerForExtent(wipeExt);

    switch (this->Wipe)
      {
      case VTK_WIPE_QUAD:
        whichInput = 0;
        break;
      case VTK_WIPE_HORIZONTAL:
        whichInput = 0;
        break;
      case VTK_WIPE_VERTICAL:
        whichInput = 0;
        break;
      case VTK_WIPE_LOWER_LEFT:
        whichInput = 0;
        break;
      case VTK_WIPE_LOWER_RIGHT:
        whichInput = 1;
        break;
      case VTK_WIPE_UPPER_LEFT:
        whichInput = 1;
        break;
      case VTK_WIPE_UPPER_RIGHT:
        whichInput = 1;
        break;
      }
    inPtr = inData[whichInput]->GetScalarPointerForExtent(wipeExt);
    switch (inData[0]->GetScalarType())
      {
      vtkTemplateMacro7(vtkImageRectilinearWipeExecute2, this,
                        inData[whichInput], (VTK_TT *)(inPtr),
                        outData, (VTK_TT *)(outPtr),
                        wipeExt, id);
      default:
        vtkErrorMacro(<< "Execute: Unknown ScalarType");
        return;
      }
    }
  // lower right
  wipeExt[0] = wholeExt[0] + this->Position[0] + 1;
  wipeExt[1] = wholeExt[1];
  wipeExt[2] = wholeExt[2];
  wipeExt[3] = wholeExt[2] + this->Position[1];
  wipeExt[4] = wholeExt[4];
  wipeExt[5] = wholeExt[5];
  if (vtkImageRectilinearWipeClampExtents(wipeExt, outExt))
    {
    switch (this->Wipe)
      {
      case VTK_WIPE_QUAD:
        whichInput = 1;
        break;
      case VTK_WIPE_HORIZONTAL:
        whichInput = 1;
        break;
      case VTK_WIPE_VERTICAL:
        whichInput = 0;
        break;
      case VTK_WIPE_LOWER_LEFT:
        whichInput = 1;
        break;
      case VTK_WIPE_LOWER_RIGHT:
        whichInput = 0;
        break;
      case VTK_WIPE_UPPER_LEFT:
        whichInput = 1;
        break;
      case VTK_WIPE_UPPER_RIGHT:
        whichInput = 1;
        break;
      }
    inPtr = inData[whichInput]->GetScalarPointerForExtent(wipeExt);
    outPtr = outData->GetScalarPointerForExtent(wipeExt);
    switch (inData[0]->GetScalarType())
      {
      vtkTemplateMacro7(vtkImageRectilinearWipeExecute2, this,
                        inData[whichInput], (VTK_TT *)(inPtr),
                        outData, (VTK_TT *)(outPtr),
                        wipeExt, id);
      default:
        vtkErrorMacro(<< "Execute: Unknown ScalarType");
        return;
      }
    }
  // upper left
  wipeExt[0] = wholeExt[0];
  wipeExt[1] = wholeExt[0] + this->Position[0];
  wipeExt[2] = wholeExt[2] + this->Position[1] + 1;
  wipeExt[3] = wholeExt[3];
  wipeExt[4] = wholeExt[4];
  wipeExt[5] = wholeExt[5];
  if (vtkImageRectilinearWipeClampExtents(wipeExt, outExt))
    {

    switch (this->Wipe)
      {
      case VTK_WIPE_QUAD:
        whichInput = 1;
        break;
      case VTK_WIPE_HORIZONTAL:
        whichInput = 0;
        break;
      case VTK_WIPE_VERTICAL:
        whichInput = 1;
        break;
      case VTK_WIPE_LOWER_LEFT:
        whichInput = 1;
        break;
      case VTK_WIPE_LOWER_RIGHT:
        whichInput = 1;
        break;
      case VTK_WIPE_UPPER_LEFT:
        whichInput = 0;
        break;
      case VTK_WIPE_UPPER_RIGHT:
        whichInput = 1;
        break;
      }
    inPtr = inData[whichInput]->GetScalarPointerForExtent(wipeExt);
    outPtr = outData->GetScalarPointerForExtent(wipeExt);
    switch (inData[0]->GetScalarType())
      {
      vtkTemplateMacro7(vtkImageRectilinearWipeExecute2, this,
                        inData[whichInput], (VTK_TT *)(inPtr),
                        outData, (VTK_TT *)(outPtr),
                        wipeExt, id);
      default:
        vtkErrorMacro(<< "Execute: Unknown ScalarType");
        return;
      }
    }
  // upper right
  wipeExt[0] = wholeExt[0] + this->Position[0] + 1;
  wipeExt[1] = wholeExt[1];
  wipeExt[2] = wholeExt[2] + this->Position[1] + 1;
  wipeExt[3] = wholeExt[3];
  wipeExt[4] = wholeExt[4];
  wipeExt[5] = wholeExt[5];
  if (vtkImageRectilinearWipeClampExtents(wipeExt, outExt))
    {
    switch (this->Wipe)
      {
      case VTK_WIPE_QUAD:
        whichInput = 0;
        break;
      case VTK_WIPE_HORIZONTAL:
        whichInput = 1;
        break;
      case VTK_WIPE_VERTICAL:
        whichInput = 1;
        break;
      case VTK_WIPE_LOWER_LEFT:
        whichInput = 1;
        break;
      case VTK_WIPE_LOWER_RIGHT:
        whichInput = 1;
        break;
      case VTK_WIPE_UPPER_LEFT:
        whichInput = 1;
        break;
      case VTK_WIPE_UPPER_RIGHT:
        whichInput = 0;
        break;
      }
    inPtr = inData[whichInput]->GetScalarPointerForExtent(wipeExt);
    outPtr = outData->GetScalarPointerForExtent(wipeExt);
    switch (inData[0]->GetScalarType())
      {
      vtkTemplateMacro7(vtkImageRectilinearWipeExecute2, this,
                        inData[whichInput], (VTK_TT *)(inPtr),
                        outData, (VTK_TT *)(outPtr),
                        wipeExt, id);
      default:
        vtkErrorMacro(<< "Execute: Unknown ScalarType");
        return;
      }
    }
}

void vtkImageRectilinearWipe::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
  os << indent << "Position: (" << this->Position[0] << ", "
     << this->Position[1] << ")\n";
  os << indent << "Wipe: ";
  switch (this->Wipe)
    {
    case VTK_WIPE_QUAD:
      os << "Quad" << endl;
      break;
    case VTK_WIPE_HORIZONTAL:
      os << "Horizontal" << endl;
      break;
    case VTK_WIPE_VERTICAL:
      os << "Vertical" << endl;
      break;
    case VTK_WIPE_LOWER_LEFT:
      os << "LowerLeft" << endl;
      break;
    case VTK_WIPE_LOWER_RIGHT:
      os << "LowerRight" << endl;
      break;
    case VTK_WIPE_UPPER_LEFT:
      os << "UpperLeft" << endl;
      break;
    case VTK_WIPE_UPPER_RIGHT:
      os << "UpperRight" << endl;
      break;
    }
}

