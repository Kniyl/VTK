/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkImageScatterPlotFilter.cc
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

Copyright (c) 1993-1995 Ken Martin, Will Schroeder, Bill Lorensen.

This software is copyrighted by Ken Martin, Will Schroeder and Bill Lorensen.
The following terms apply to all files associated with the software unless
explicitly disclaimed in individual files. This copyright specifically does
not apply to the related textbook "The Visualization Toolkit" ISBN
013199837-4 published by Prentice Hall which is covered by its own copyright.

The authors hereby grant permission to use, copy, and distribute this
software and its documentation for any purpose, provided that existing
copyright notices are retained in all copies and that this notice is included
verbatim in any distributions. Additionally, the authors grant permission to
modify this software and its documentation for any purpose, provided that
such modifications are not distributed without the explicit consent of the
authors and that existing copyright notices are retained in all copies. Some
of the algorithms implemented by this software are patented, observe all
applicable patent law.

IN NO EVENT SHALL THE AUTHORS OR DISTRIBUTORS BE LIABLE TO ANY PARTY FOR
DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES ARISING OUT
OF THE USE OF THIS SOFTWARE, ITS DOCUMENTATION, OR ANY DERIVATIVES THEREOF,
EVEN IF THE AUTHORS HAVE BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

THE AUTHORS AND DISTRIBUTORS SPECIFICALLY DISCLAIM ANY WARRANTIES, INCLUDING,
BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
PARTICULAR PURPOSE, AND NON-INFRINGEMENT.  THIS SOFTWARE IS PROVIDED ON AN
"AS IS" BASIS, AND THE AUTHORS AND DISTRIBUTORS HAVE NO OBLIGATION TO PROVIDE
MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.


=========================================================================*/
#include "vtkImageScatterPlotFilter.hh"



//----------------------------------------------------------------------------
vtkImageScatterPlotFilter::vtkImageScatterPlotFilter()
{
  this->SetAxes4d(VTK_IMAGE_X_AXIS, VTK_IMAGE_Y_AXIS, 
		  VTK_IMAGE_Z_AXIS, VTK_IMAGE_COMPONENT_AXIS);
  this->SetOutputDataType(VTK_IMAGE_UNSIGNED_SHORT);
  // set up for a 2d (256x256) image representing 0->256 in the components.
  this->ImageRegion.SetBounds2d(0, 255, 0, 255);
  this->AspectRatio = 1.0;
  // We want to request our input all by ourself.
  this->UseExecuteMethodOff();
  this->InRegion.SetBounds4d(0, 255, 0, 255, 47, 47, 0, 1);
  
}


//----------------------------------------------------------------------------
// Description:
// Set the input of the filter, sets the default bounds of the InRegion as
// the ImageBounds of the input.
void vtkImageScatterPlotFilter::SetInput(vtkImageSource *input)
{
  this->vtkImageFilter::SetInput(input);
  input->UpdateImageInformation(&(this->InRegion));
}


//----------------------------------------------------------------------------
// Description:
// Set coordinate system of the filter.
void vtkImageScatterPlotFilter::SetAxes(int *axes)
{
  this->vtkImageCachedSource::SetAxes(axes);
  this->InRegion.SetAxes(axes);
  this->ImageRegion.SetAxes(axes);
}


//----------------------------------------------------------------------------
// Description:
// Just sets the ImageBounds to be the bounds of the OutRegion.
void vtkImageScatterPlotFilter::ComputeOutputImageInformation(
		    vtkImageRegion *inRegion, vtkImageRegion *outRegion)
{
  // Avoid warnings
  inRegion = inRegion;
  outRegion->SetImageBounds(this->ImageRegion.GetBounds());
}


  
//----------------------------------------------------------------------------
// Description:
// The Required input bounds are the bounds of the InRegion.
// Note used in this implementation
void vtkImageScatterPlotFilter::ComputeRequiredInputRegionBounds(
				       vtkImageRegion *outRegion, 
				       vtkImageRegion *inRegion)
{
  outRegion = outRegion;
  inRegion->SetBounds(this->InRegion.GetBounds());
}


//----------------------------------------------------------------------------
// Description:
// This templated function executes the filter for any type of data.
template <class T>
void vtkImageScatterPlotFilterUpdate(vtkImageScatterPlotFilter *self,
			     vtkImageRegion *inRegion, T *inPtr,
			     vtkImageRegion *outRegion, unsigned short *outPtr)
{
  int outMin0, outMax0, outMin1, outMax1;
  int inMin0, inMax0, inMin1, inMax1;
  int inInc0, inInc1, inInc2, inInc3;
  int idx0, idx1;
  T *inPtr0, *inPtr1;
  float aspectRatio = self->GetAspectRatio();
  int coordinates[2];
  
  inRegion->GetBounds2d(inMin0, inMax0, inMin1, inMax1);
  inRegion->GetIncrements4d(inInc0, inInc1, inInc2, inInc3);
  outRegion->GetBounds2d(outMin0, outMax0, outMin1, outMax1);
  
  // loop over all the input pixels (2d images).
  inPtr1 = inPtr;
  for (idx1 = inMin1; idx1 <= inMax1; ++idx1)
    {
    inPtr0 = inPtr1;
    for (idx0 = inMin0; idx0 <= inMax0; ++idx0)
      {
      // get and convert components into indicies
      coordinates[0] = (int)((float)(*inPtr0) / aspectRatio);
      coordinates[1] = (int)((float)(*(inPtr0+inInc3)) / aspectRatio);
      // Add to scatter plot
      if (outMin0 <= coordinates[0] && coordinates[0] <= outMax0 &&
	  outMin1 <= coordinates[1] && coordinates[1] <= outMax1)
	{
	outPtr = (unsigned short *)(outRegion->GetVoidPointer2d(coordinates));
	++(*outPtr);
	}

      inPtr0 += inInc0;
      }
    inPtr1 += inInc1;
    }
}



//----------------------------------------------------------------------------
// Description:
// This function requests the input as 2d images. It does not handle
// failed input requests at all.
void vtkImageScatterPlotFilter::UpdateRegion(vtkImageRegion *outRegion)
{
  // only handle 4d.
  vtkImageRegion *inRegion = NULL;
  void *inPtr, *outPtr;
  int bounds[8];
  int min0, max0, min1, max1, min2, max2;
  int inc0, inc1;
  int idx0, idx1, idx2;
  unsigned short *outPtr0, *outPtr1;

  
  // We need an input (error checking)
  if ( ! this->Input)
    {
    vtkErrorMacro(<< "UpdateRegion4d: Input not set.");
    }
  
  // this filter expects that output is unsigned short
  if (outRegion->GetDataType() != VTK_IMAGE_UNSIGNED_SHORT)
    {
    vtkErrorMacro(<< "Execute: output must be unsigned short");
    return;
    }

  // Allocate output region (input is requested in pieces.)
  this->Output->AllocateRegion(outRegion);

  // Set all the output pixels to 0;
  outRegion->GetBounds2d(min0, max0, min1, max1);
  outRegion->GetIncrements2d(inc0, inc1);
  outPtr1 = (unsigned short *)(outRegion->GetVoidPointer2d());
  for (idx1 = min1; idx1 <= max1; ++idx1)
    {
    outPtr0 = outPtr1;
    for (idx0 = min0; idx0 <= max0; ++idx0)
      {
      *outPtr0 = 0;
      
      outPtr0 += inc0;
      }
    outPtr1 += inc1;
    }
  
  // Input region can only have two components for now.
  this->InRegion.GetBounds4d(bounds);
  if (bounds[7] - bounds[6] != 1)
    {
    vtkErrorMacro(<< "Update Region: I only generate 2d plots.");
    return;
    }
  
  // loop over the multi spectral input images. (third axis)
  min2 = bounds[4];
  max2 = bounds[5];
  for (idx2 = min2; idx2 <= max2; ++idx2)
    {
    // Set up the region to have the bounds of one 2d image.
    bounds[4] = bounds[5] = idx2;
    // Request Image.
    inRegion = this->GetInputRegion4d(bounds);
    if ( ! inRegion->IsAllocated())
      {
      vtkErrorMacro(<< "UpdateRegion:Input request failed.");
      return;
      }
    // Call the template function to add to this scatter plot.
    inPtr = inRegion->GetVoidPointer4d();
    outPtr = inRegion->GetVoidPointer4d();
    switch (inRegion->GetDataType())
      {
      case VTK_IMAGE_FLOAT:
	vtkImageScatterPlotFilterUpdate(this, inRegion, (float *)(inPtr), 
					outRegion, (unsigned short *)(outPtr));
	break;
      case VTK_IMAGE_INT:
	vtkImageScatterPlotFilterUpdate(this, inRegion, (int *)(inPtr), 
					outRegion, (unsigned short *)(outPtr));
	break;
      case VTK_IMAGE_SHORT:
	vtkImageScatterPlotFilterUpdate(this, inRegion, (short *)(inPtr), 
					outRegion, (unsigned short *)(outPtr));
	break;
      case VTK_IMAGE_UNSIGNED_SHORT:
	vtkImageScatterPlotFilterUpdate(this, 
					inRegion, (unsigned short *)(inPtr),
					outRegion, (unsigned short *)(outPtr));
	break;
      case VTK_IMAGE_UNSIGNED_CHAR:
	vtkImageScatterPlotFilterUpdate(this, 
					inRegion, (unsigned char *)(inPtr), 
					outRegion, (unsigned short *)(outPtr));
	break;
      default:
	vtkErrorMacro(<< "Execute: Unknown DataType");
	return;
      }
    // We are done with this input, so delete it.
    inRegion->Delete();
    inRegion = NULL;
    }
}

  














