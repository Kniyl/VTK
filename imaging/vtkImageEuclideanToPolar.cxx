/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkImageEuclideanToPolar.cxx
  Language:  C++
  Date:      $Date$
  Version:   $Revision$
  Thanks:    Thanks to C. Charles Law who developed this class.

Copyright (c) 1993-2000 Ken Martin, Will Schroeder, Bill Lorensen.

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
#include "vtkImageEuclideanToPolar.h"
#include "vtkObjectFactory.h"



//------------------------------------------------------------------------------
vtkImageEuclideanToPolar* vtkImageEuclideanToPolar::New()
{
  // First try to create the object from the vtkObjectFactory
  vtkObject* ret = vtkObjectFactory::CreateInstance("vtkImageEuclideanToPolar");
  if(ret)
    {
    return (vtkImageEuclideanToPolar*)ret;
    }
  // If the factory was unable to create the object, then create it here.
  return new vtkImageEuclideanToPolar;
}






//----------------------------------------------------------------------------
vtkImageEuclideanToPolar::vtkImageEuclideanToPolar()
{
  this->ThetaMaximum = 255.0;
}

//----------------------------------------------------------------------------
// This templated function executes the filter for any type of data.
template <class T>
static void vtkImageEuclideanToPolarExecute(vtkImageEuclideanToPolar *self,
				    vtkImageData *inData, T *inPtr,
				    vtkImageData *outData, T *outPtr,
				    int outExt[6], int id)
{
  int idxX, idxY, idxZ;
  int maxC, maxX, maxY, maxZ;
  int inIncX, inIncY, inIncZ;
  int outIncX, outIncY, outIncZ;
  unsigned long count = 0;
  unsigned long target;
  float X, Y, Theta, R;
  float thetaMax = self->GetThetaMaximum();
  
  // find the region to loop over
  maxC = inData->GetNumberOfScalarComponents();
  maxX = outExt[1] - outExt[0]; 
  maxY = outExt[3] - outExt[2]; 
  maxZ = outExt[5] - outExt[4];
  target = (unsigned long)((maxZ+1)*(maxY+1)/50.0);
  target++;
  
  // Get increments to march through data 
  inData->GetContinuousIncrements(outExt, inIncX, inIncY, inIncZ);
  outData->GetContinuousIncrements(outExt, outIncX, outIncY, outIncZ);

  // Loop through ouput pixels
  for (idxZ = 0; idxZ <= maxZ; idxZ++)
    {
    for (idxY = 0; !self->AbortExecute && idxY <= maxY; idxY++)
      {
      if (!id) 
	{
	if (!(count%target))
	  {
	  self->UpdateProgress(count/(50.0*target));
	  }
	count++;
	}
      for (idxX = 0; idxX <= maxX; idxX++)
	{
	// Pixel operation
	X = (float)(*inPtr);
	Y = (float)(inPtr[1]);

	if ((X == 0.0) && (Y == 0.0))
	  {
	  Theta = 0.0;
	  R = 0.0;
	  }
	else
	  {
	  Theta = atan2(Y, X) * thetaMax / 6.2831853;
	  if (Theta < 0.0)
	    {
	    Theta += thetaMax;
	    }
	  R = sqrt(X*X + Y*Y);
	  }
	
	*outPtr = (T)(Theta);
	outPtr[1] = (T)(R);
	inPtr += maxC;
	outPtr += maxC;
	}
      outPtr += outIncY;
      inPtr += inIncY;
      }
    outPtr += outIncZ;
    inPtr += inIncZ;
    }
}

//----------------------------------------------------------------------------
void vtkImageEuclideanToPolar::ThreadedExecute(vtkImageData *inData, 
				       vtkImageData *outData,
				       int outExt[6], int id)
{
  void *inPtr = inData->GetScalarPointerForExtent(outExt);
  void *outPtr = outData->GetScalarPointerForExtent(outExt);
  
  vtkDebugMacro(<< "Execute: inData = " << inData 
		<< ", outData = " << outData);
  
  // this filter expects that input is the same type as output.
  if (inData->GetScalarType() != outData->GetScalarType())
    {
    vtkErrorMacro(<< "Execute: input ScalarType, " << inData->GetScalarType()
    << ", must match out ScalarType " << outData->GetScalarType());
    return;
    }
  
  // input must have at least two components
  if (inData->GetNumberOfScalarComponents() < 2)
    {
    vtkErrorMacro(<< "Execute: input does not have at least two components");
    return;
    }

  switch (inData->GetScalarType())
    {
    case VTK_DOUBLE:
      vtkImageEuclideanToPolarExecute(this, 
			      inData, (double *)(inPtr), 
			      outData, (double *)(outPtr), outExt, id);
      break;
    case VTK_FLOAT:
      vtkImageEuclideanToPolarExecute(this, 
			      inData, (float *)(inPtr), 
			      outData, (float *)(outPtr), outExt, id);
      break;
    case VTK_LONG:
      vtkImageEuclideanToPolarExecute(this, 
			      inData, (long *)(inPtr), 
			      outData, (long *)(outPtr), outExt, id);
      break;
    case VTK_UNSIGNED_LONG:
      vtkImageEuclideanToPolarExecute(this, 
			      inData, (unsigned long *)(inPtr), 
			      outData, (unsigned long *)(outPtr), 
			      outExt, id);
      break;
    case VTK_INT:
      vtkImageEuclideanToPolarExecute(this, 
			      inData, (int *)(inPtr), 
			      outData, (int *)(outPtr), outExt, id);
      break;
    case VTK_UNSIGNED_INT:
      vtkImageEuclideanToPolarExecute(this, 
			      inData, (unsigned int *)(inPtr), 
			      outData, (unsigned int *)(outPtr), 
			      outExt, id);
      break;
    case VTK_SHORT:
      vtkImageEuclideanToPolarExecute(this, 
			      inData, (short *)(inPtr), 
			      outData, (short *)(outPtr), outExt, id);
      break;
    case VTK_UNSIGNED_SHORT:
      vtkImageEuclideanToPolarExecute(this, 
			      inData, (unsigned short *)(inPtr), 
			      outData, (unsigned short *)(outPtr), 
			      outExt, id);
      break;
    case VTK_CHAR:
      vtkImageEuclideanToPolarExecute(this, 
			      inData, (char *)(inPtr), 
			      outData, (char *)(outPtr), outExt, id);
      break;
    case VTK_UNSIGNED_CHAR:
      vtkImageEuclideanToPolarExecute(this, 
			      inData, (unsigned char *)(inPtr), 
			      outData, (unsigned char *)(outPtr), 
			      outExt, id);
      break;
    default:
      vtkErrorMacro(<< "Execute: Unknown ScalarType");
      return;
    }
}

void vtkImageEuclideanToPolar::PrintSelf(ostream& os, vtkIndent indent)
{
  vtkImageToImageFilter::PrintSelf(os,indent);

  os << indent << "Maximum Angle: " << this->ThetaMaximum << "\n";
}
