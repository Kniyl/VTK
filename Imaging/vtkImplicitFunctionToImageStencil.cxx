/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkImplicitFunctionToImageStencil.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkImplicitFunctionToImageStencil.h"

#include "vtkImageStencilData.h"
#include "vtkImplicitFunction.h"
#include "vtkObjectFactory.h"

#include <math.h>

vtkCxxRevisionMacro(vtkImplicitFunctionToImageStencil, "1.9");
vtkStandardNewMacro(vtkImplicitFunctionToImageStencil);
vtkCxxSetObjectMacro(vtkImplicitFunctionToImageStencil,Input, vtkImplicitFunction);

//----------------------------------------------------------------------------
vtkImplicitFunctionToImageStencil::vtkImplicitFunctionToImageStencil()
{
  this->Threshold = 0;
  this->Input = NULL;
}

//----------------------------------------------------------------------------
vtkImplicitFunctionToImageStencil::~vtkImplicitFunctionToImageStencil()
{
  this->SetInput(NULL);
}

//----------------------------------------------------------------------------
void vtkImplicitFunctionToImageStencil::PrintSelf(ostream& os,
                                                  vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "Input: " << this->Input << "\n";
  os << indent << "Threshold: " << this->Threshold << "\n";
}

//----------------------------------------------------------------------------
// set up the clipping extents from an implicit function by brute force
// (i.e. by evaluating the function at each and every voxel)
void vtkImplicitFunctionToImageStencil::ThreadedExecute(vtkImageStencilData 
                                                        *data,
                                                        int extent[6], int id)
{
  vtkImplicitFunction *function = this->Input;
  double *spacing = data->GetSpacing();
  double *origin = data->GetOrigin();
  double threshold = this->Threshold;

  // if the input is not set then punt
  if (!function)
    {
    return;
    }

  // for conversion of (idX,idY,idZ) into (x,y,z)
  double point[3];

  // for keeping track of progress
  unsigned long count = 0;
  unsigned long target = (unsigned long)
    ((extent[5] - extent[4] + 1)*(extent[3] - extent[2] + 1)/50.0);
  target++;

  // loop through all voxels
  for (int idZ = extent[4]; idZ <= extent[5]; idZ++)
    {
    point[2] = idZ*spacing[2] + origin[2];

    for (int idY = extent[2]; idY <= extent[3]; idY++)
      {
      point[1] = idY*spacing[1] + origin[1];
      int state = 1; // inside or outside, start outside
      int r1 = extent[0];
      int r2 = extent[1];

      if (id == 0)
        { // update progress if we're the main thread
        if (count%target == 0) 
          {
          this->UpdateProgress(count/(50.0*target));
          }
        count++;
        }

      for (int idX = extent[0]; idX <= extent[1]; idX++)
        {
        point[0] = idX*spacing[0] + origin[0];
        int newstate = 1;
        if (function->FunctionValue(point) < threshold)
          {
          newstate = -1;
          if (newstate != state)
            { // sub extent starts
            r1 = idX;
            }
          }
        else if (newstate != state)
          { // sub extent ends
          r2 = idX - 1;
          data->InsertNextExtent(r1, r2, idY, idZ);
          }
        state = newstate;
        } // for idX
      if (state == -1)
        { // if inside at end, cap off the sub extent
        data->InsertNextExtent(r1, extent[1], idY, idZ);
        }
      } // for idY    
    } // for idZ
}

void vtkImplicitFunctionToImageStencil::ExecuteInformation()
{
  vtkImageStencilData *output;
  
  output = this->GetOutput();
  
  // this is an odd source that can produce any requested size.  so its whole
  // extent is essentially infinite. This would not be a great source to
  // connect to some sort of writer or viewer. For a sanity check we will
  // limit the size produced to something reasonable (depending on your
  // definition of reasonable)
  output->SetWholeExtent(0, VTK_LARGE_INTEGER >> 2,
                         0, VTK_LARGE_INTEGER >> 2,
                         0, VTK_LARGE_INTEGER >> 2);
}
