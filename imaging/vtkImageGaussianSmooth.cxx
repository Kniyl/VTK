/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkImageGaussianSmooth.cxx
  Language:  C++
  Date:      $Date$
  Version:   $Revision$
  Thanks:    Thanks to C. Charles Law who developed this class.

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
#include <math.h>
#include "vtkImageGaussianSmooth.h"

//----------------------------------------------------------------------------
vtkImageGaussianSmooth::vtkImageGaussianSmooth()
{
  this->Radius = -1;
  this->Kernel = this->TempKernel = NULL;
  this->Dimensionality = 1; // note: this overrides Standard deviation.
  this->StandardDeviations[0] = 2.0;
  this->StandardDeviations[1] = 2.0;
  this->StandardDeviations[2] = 2.0;
  this->RadiusFactors[0] = 1.5;
  this->RadiusFactors[1] = 1.5;
  this->RadiusFactors[2] = 1.5;
}

//----------------------------------------------------------------------------
vtkImageGaussianSmooth::~vtkImageGaussianSmooth()
{
  if (this->Kernel)
    {
    delete [] this->Kernel;
    this->Kernel = NULL;
    delete [] this->TempKernel;
    this->TempKernel = NULL;
    }
}

//----------------------------------------------------------------------------
void vtkImageGaussianSmooth::PrintSelf(ostream& os, vtkIndent indent)
{
  // int idx;
  
  this->vtkImageFilter::PrintSelf(os, indent);
  //os << indent << "BoundaryRescale: " << this->BoundaryRescale << "\n";
}

//----------------------------------------------------------------------------
void vtkImageGaussianSmooth::ComputeKernel(float std, float factor)
{
  int idx, size;
  float sum, x;
  
  this->Radius = (int)(std * factor);
  if (this->Kernel)
    {
    delete [] this->Kernel;
    delete [] this->TempKernel;
    }

  // allocate space for kernels 
  size = 2 * this->Radius + 1;
  this->Kernel = new float[size];
  this->TempKernel = new float[size];
  
  // handle special case
  if (std == 0.0)
    {
    this->Kernel[0] = 1.0;
    return;
    }
  
  // fill in kernel
  sum = 0.0;
  for (idx = 0; idx < size; ++idx)
    {
    x = idx - this->Radius;
    sum += this->Kernel[idx] = 
      exp(- (x*x) / (2.0 * std * std));
    }

  // normalize
  sum = 1.0 / sum;
  for (idx = 0; idx < size; ++idx)
    {
    this->Kernel[idx] *= sum;
    }
}
  

//----------------------------------------------------------------------------
void vtkImageGaussianSmooth::ExecuteImageInformation()
{
  this->Output->SetScalarType(VTK_FLOAT);
}

//----------------------------------------------------------------------------
void vtkImageGaussianSmooth::ComputeRequiredInputUpdateExtent(int inExt[6], 
							      int outExt[6])
{
  int *wholeExtent;
  int idx, radius;

  // copy
  memcpy((void *)inExt, (void *)outExt, 6 * sizeof(int));
  // Expand filtered axes
  wholeExtent = this->Input->GetWholeExtent();
  for (idx = 0; idx < this->Dimensionality; ++idx)
    {
    radius = (int)(this->StandardDeviations[idx] * this->RadiusFactors[idx]);
    inExt[idx*2] -= radius;
    if (inExt[idx*2] < wholeExtent[idx*2])
      {
      inExt[idx*2] = wholeExtent[idx*2];
      }

    inExt[idx*2+1] += radius;
    if (inExt[idx*2+1] > wholeExtent[idx*2+1])
      {
      inExt[idx*2+1] = wholeExtent[idx*2+1];
      }
    }
}


//----------------------------------------------------------------------------
// For a given position along the convolution axis, this method loops over 
// all other axes, and performs the convolution. Boundary conditions handled
// previously.
template <class T>
static void 
vtkImageGaussianSmoothExecute(vtkImageGaussianSmooth *self, int axis,
		      float *kernel, int kernelSize,
		      vtkImageData *inData, T *inPtrC,
		      vtkImageData *outData, int outExt[6], float *outPtrC)
{
  int maxC, max0, max1;
  int idxC, idx0, idx1, idxK;
  int *inIncs, *outIncs;
  int inInc0, inInc1, inIncK, outInc0, outInc1;
  float *outPtr1, *outPtr0;
  T *inPtr1, *inPtr0, *inPtrK;
  float *ptrK, sum;
  
  // avoid warnings
  self = self;
  max0 = max1 = inInc0 = inInc1 = outInc0 = outInc1 = 0;
  
  // I am counting on the fact that tight loops (component on outside)
  // is more important than cache misses from shuffled access.

  // Do the correct shuffling of the axes (increments, extents)
  inIncs = inData->GetIncrements();
  outIncs = outData->GetIncrements();
  inIncK = inIncs[axis];
  maxC = outData->GetNumberOfScalarComponents();
  switch (axis)
    {
    case 0:
      inInc0 = inIncs[1];  inInc1 = inIncs[2];
      outInc0 = outIncs[1];  outInc1 = outIncs[2];
      max0 = outExt[3] - outExt[2] + 1;   max1 = outExt[5] - outExt[4] + 1;
      break;
    case 1:
      inInc0 = inIncs[0];  inInc1 = inIncs[2];
      outInc0 = outIncs[0];  outInc1 = outIncs[2];
      max0 = outExt[1] - outExt[0] + 1;   max1 = outExt[5] - outExt[4] + 1;
      break;
    case 2:
      inInc0 = inIncs[0];  inInc1 = inIncs[1];
      outInc0 = outIncs[0];  outInc1 = outIncs[1];
      max0 = outExt[1] - outExt[0] + 1;   max1 = outExt[3] - outExt[2] + 1;
      break;
    }
  
  for (idxC = 0; idxC < maxC; ++idxC)
    {
    inPtr1 = inPtrC;
    outPtr1 = outPtrC;    
    for (idx1 = 0; idx1 < max1; ++idx1)
      {
      inPtr0 = inPtr1;
      outPtr0 = outPtr1;    
      for (idx0 = 0; idx0 < max0; ++idx0)
	{
	inPtrK = inPtr0;
	ptrK = kernel;
	sum = 0.0;
	// too bad this short loop has to be the inner most loop
	for (idxK = 0; idxK < kernelSize; ++idxK)
	  {
	  sum += *ptrK * (float)(*inPtrK);
	  ++ptrK;
	  inPtrK += inIncK;
	  }
	*outPtr0 = sum;
	inPtr0 += inInc0;
	outPtr0 += outInc0;
	}
      inPtr1 += inInc1;
      outPtr1 += outInc1;
      }
    ++inPtrC;
    ++outPtrC;
    }
}


//----------------------------------------------------------------------------
// Description:
// This method convolves over one axis. It loops over the convolved axis,
// and handles boundary conditions.
void vtkImageGaussianSmooth::ExecuteAxis(int axis,
					 vtkImageData *inData, int inExt[6],
					 vtkImageData *outData, int outExt[6])
{
  int idxA, idx, max;
  int *wholeExtent, wholeMax, wholeMin;
  float *kernel, sum;
  int kernelSize, kernelLeftClip, kernelRightClip;
  void *inPtr;
  float *outPtr;
  int coords[3], *outIncs, outIncA;

  // Compute the correct kernel for this axis (also set Radius)
  this->ComputeKernel(this->StandardDeviations[axis], 
		      this->RadiusFactors[axis]);
  
  // Get the correct starting pointer of the output
  outPtr = (float *)outData->GetScalarPointerForExtent(outExt);
  outIncs = outData->GetIncrements();
  outIncA = outIncs[axis];

  // Determine default starting position of input
  coords[0] = inExt[0];
  coords[1] = inExt[2];
  coords[2] = inExt[4];
  
  // get whole extent for boundary checking ...
  wholeExtent = this->Input->GetWholeExtent();
  wholeMin = wholeExtent[axis*2];
  wholeMax = wholeExtent[axis*2+1];  

  // loop over the convolution axis
  max = outExt[axis*2+1];
  for (idxA = outExt[axis*2]; idxA <= max; ++idxA, outPtr += outIncA)
    {
    // left boundary condition
    coords[axis] = idxA - this->Radius;
    kernelLeftClip = wholeMin - coords[axis];
    if (kernelLeftClip > 0)
      { // front of kernel is cut off ("kernelStart" samples)
      coords[axis] += kernelLeftClip;
      }
    else
      {
      kernelLeftClip = 0;
      }
    // Right boundary condition
    kernelRightClip = (idxA + this->Radius) - wholeMax;
    if (kernelRightClip < 0)
      {
      kernelRightClip = 0;
      }
    
    // compute the clipped Kernel.
    if (kernelLeftClip == 0 && kernelRightClip == 0)
      { // kernel not clipped
      kernel = this->Kernel;
      kernelSize = 2 * this->Radius + 1;
      }
    else
      {
      // copy clipped kernel
      sum = 0.0;
      kernelSize = 2 * this->Radius + 1 - kernelRightClip - kernelLeftClip;
      for (idx = 0; idx < kernelSize; ++idx)
	{
	sum += this->TempKernel[idx] = this->Kernel[idx + kernelLeftClip];
	}
      // normalize
      sum = 1.0 / sum;
      for (idx = 0; idx < kernelSize; ++idx)
	{
	this->TempKernel[idx] *= sum;
	}
      kernel = this->TempKernel;
      }
    
    /* now do the convolution on the rest of the axes */
    inPtr = inData->GetScalarPointer(coords);
    switch (inData->GetScalarType())
      {
      case VTK_FLOAT:
	vtkImageGaussianSmoothExecute(this, axis, kernel, kernelSize,
				      inData, (float *)(inPtr),
				      outData, outExt, outPtr);
	break;
      case VTK_INT:
	vtkImageGaussianSmoothExecute(this, axis, kernel, kernelSize,
				      inData, (int *)(inPtr),
				      outData, outExt, outPtr);
	break;
      case VTK_SHORT:
	vtkImageGaussianSmoothExecute(this, axis, kernel, kernelSize,
				      inData, (short *)(inPtr),
				      outData, outExt, outPtr);
	break;
      case VTK_UNSIGNED_SHORT:
	vtkImageGaussianSmoothExecute(this, axis, kernel, kernelSize,
				      inData, (unsigned short *)(inPtr),
				      outData, outExt, outPtr);
	break;
      case VTK_UNSIGNED_CHAR:
	vtkImageGaussianSmoothExecute(this, axis, kernel, kernelSize,
				      inData, (unsigned char *)(inPtr),
				      outData, outExt, outPtr);
	break;
      default:
	vtkErrorMacro("Unknown scalar type");
	return;
      }
    }
}



  
//----------------------------------------------------------------------------
// Description:
// This method decomposes the gaussian and smooths along each axis.
void vtkImageGaussianSmooth::ThreadedExecute(vtkImageData *inData, 
					     vtkImageData *outData,
					     int outExt[6], int id)
{
  int inExt[6];
  
  id = id;
  
  // this filter expects the output to be type float
  if (outData->GetScalarType() != VTK_FLOAT)
    {
    vtkErrorMacro("Execute: output ScalarType, " 
		  << outData->GetScalarType() << ", must be float");
    return;
    }

  // Decompose
  this->ComputeRequiredInputUpdateExtent(inExt, outExt);
  switch (this->Dimensionality)
    {
    case 1:
      ExecuteAxis(0, inData, inExt, outData, outExt);
      break;
    case 2:
      int tempExt[6];
      vtkImageData *tempData;
      // compute intermediate extent
      tempExt[0] = inExt[0];  tempExt[1] = inExt[1];
      tempExt[2] = outExt[2];  tempExt[3] = outExt[3];
      tempExt[4] = inExt[4];  tempExt[5] = inExt[5];
      // create a temp data for intermediate results
      tempData = vtkImageData::New();
      tempData->SetExtent(tempExt);
      tempData->SetNumberOfScalarComponents(inData->GetNumberOfScalarComponents());
      tempData->SetScalarType(VTK_FLOAT);
      ExecuteAxis(1, inData, inExt, tempData, tempExt);
      ExecuteAxis(0, tempData, tempExt, outData, outExt);
      // release temporary data
      tempData->Delete();
      tempData = NULL;
      break;
    case 3:
      // we do z first because it is most likely smallest
      int temp0Ext[6], temp1Ext[6];
      vtkImageData *temp0Data, *temp1Data;
      // compute intermediate extents
      temp0Ext[0] = inExt[0];  temp0Ext[1] = inExt[1];
      temp0Ext[2] = inExt[2];  temp0Ext[3] = inExt[3];
      temp0Ext[4] = outExt[4];  temp0Ext[5] = outExt[5];

      temp1Ext[0] = inExt[0];  temp1Ext[1] = inExt[1];
      temp1Ext[2] = outExt[2];  temp1Ext[3] = outExt[3];
      temp1Ext[4] = outExt[4];  temp1Ext[5] = outExt[5];
      
      // create a temp data for intermediate results
      temp0Data = vtkImageData::New();
      temp0Data->SetExtent(temp0Ext);
      temp0Data->SetNumberOfScalarComponents(inData->GetNumberOfScalarComponents());
      temp0Data->SetScalarType(VTK_FLOAT);

      temp1Data = vtkImageData::New();
      temp1Data->SetExtent(temp1Ext);
      temp1Data->SetNumberOfScalarComponents(inData->GetNumberOfScalarComponents());
      temp1Data->SetScalarType(VTK_FLOAT);
      ExecuteAxis(2, inData, inExt, temp0Data, temp0Ext);
      ExecuteAxis(1, temp0Data, temp0Ext, temp1Data, temp1Ext);
      temp0Data->Delete();
      temp0Data = NULL;
      ExecuteAxis(0, temp1Data, temp1Ext, outData, outExt);
      temp1Data->Delete();
      temp1Data = NULL;
      break;
    }  
}















