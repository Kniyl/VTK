/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkImageReslice.cxx
  Language:  C++
  Date:      $Date$
  Version:   $Revision$
  Thanks:    Thanks to David G Gobbi who developed this class.

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
#include <limits.h>
#include <float.h>
#include <math.h>
#include "vtkImageCache.h"
#include "vtkImageReslice.h"

//----------------------------------------------------------------------------
vtkImageReslice::vtkImageReslice()
{
  this->OutputSpacing[0] = 1;
  this->OutputSpacing[1] = 1;
  this->OutputSpacing[2] = 1;

  this->OutputOrigin[0] = FLT_MAX; // flag to set defaults later
  this->OutputExtent[0] = INT_MAX; // ditto

  this->Interpolate = 0; // nearest-neighbor interpolation by default
  this->Optimization = 1; // optimizations seem to finally be stable...
  this->BackgroundLevel = 0;

  this->ResliceTransform = NULL;
}

//----------------------------------------------------------------------------
vtkImageReslice::~vtkImageReslice()
{
  this->SetResliceTransform(NULL);
}

//----------------------------------------------------------------------------
void vtkImageReslice::PrintSelf(ostream& os, vtkIndent indent)
{
  vtkImageFilter::PrintSelf(os,indent);

  if (this->ResliceTransform)
    {
    os << indent << "ResliceTransform: " << this->ResliceTransform << "\n";
    this->ResliceTransform->PrintSelf(os,indent.GetNextIndent());
    }
  os << indent << "OutputSpacing: " << this->OutputSpacing[0] << " " <<
    this->OutputSpacing[1] << " " << this->OutputSpacing[2] << "\n";
  os << indent << "OutputOrigin: " << this->OutputOrigin[0] << " " <<
    this->OutputOrigin[1] << " " << this->OutputOrigin[2] << "\n";
  os << indent << "OutputExtent: " << this->OutputExtent[0] << " " <<
    this->OutputExtent[1] << " " << this->OutputExtent[2] << " " <<
    this->OutputExtent[3] << " " << this->OutputExtent[4] << " " <<
    this->OutputExtent[5] << "\n";
  os << indent << "Interpolate: " << this->Interpolate << "\n";
  os << indent << "BackgroundLevel: " << this->BackgroundLevel << "\n";
}

//----------------------------------------------------------------------------
void vtkImageReslice::ComputeIndexMatrix(vtkMatrix4x4 *matrix)
{
  int i;
  float inOrigin[3];
  float inSpacing[3];
  float outOrigin[3];
  float outSpacing[3];

  this->GetInput()->GetSpacing(inSpacing);
  this->GetInput()->GetOrigin(inOrigin);
  this->GetOutput()->GetSpacing(outSpacing);
  this->GetOutput()->GetOrigin(outOrigin);  
  
  vtkTransform *transform = vtkTransform::New();
  vtkMatrix4x4 *inMatrix = vtkMatrix4x4::New();
  vtkMatrix4x4 *outMatrix = vtkMatrix4x4::New();

  if (this->GetResliceTransform())
    {
    transform->SetMatrix(*(this->GetResliceTransform()->GetMatrixPointer()));
    }
  
  // the outMatrix takes OutputData indices to OutputData coordinates,
  // the inMatrix takes InputData coordinates to InputData indices

  for (i = 0; i < 3; i++) 
    {
    inMatrix->Element[i][i] = 1/inSpacing[i];
    inMatrix->Element[i][3] = -inOrigin[i]/inSpacing[i];
    outMatrix->Element[i][i] = outSpacing[i];
    outMatrix->Element[i][3] = outOrigin[i];
    }
  
  transform->PreMultiply();
  transform->Concatenate(outMatrix);
  transform->PostMultiply();
  transform->Concatenate(inMatrix);

  transform->GetMatrix(matrix);
  
  transform->Delete();
  inMatrix->Delete();
  outMatrix->Delete();
}

//----------------------------------------------------------------------------
void vtkImageReslice::ComputeRequiredInputUpdateExtent(int inExt[6], 
						       int outExt[6])
{
  int i,j,k;
  int idX,idY,idZ;
  vtkMatrix4x4 *matrix = vtkMatrix4x4::New();
  float *xAxis, *yAxis, *zAxis, *origin;
  double point[4],w;

  // break transform down (to match vtkImageResliceExecute)
  this->ComputeIndexMatrix(matrix);
  matrix->Transpose();
  xAxis = matrix->Element[0];
  yAxis = matrix->Element[1];
  zAxis = matrix->Element[2];
  origin = matrix->Element[3];

  for (i = 0; i < 3; i++)
    {
    inExt[2*i] = INT_MAX;
    inExt[2*i+1] = INT_MIN;
    }
  
  for (i = 0; i < 8; i++)
    {
    // calculate transform using method in vtkImageResliceExecute
    idX = outExt[i%2];
    idY = outExt[2+(i/2)%2];
    idZ = outExt[4+(i/4)%2];
    
    for (j = 0; j < 4; j++) 
      {
      point[j] = origin[j] + idZ*zAxis[j];
      point[j] = point[j] + idY*yAxis[j];
      }
    
    w = point[3] + idX*xAxis[3];
    
    if (this->Interpolate)
      {
      for (j = 0; j < 3; j++) 
	{
	k = int(float((point[j]+idX*xAxis[j])/w));
	if (k < inExt[2*j]) inExt[2*j] = k; 
	k = int(float((point[j]+idX*xAxis[j])/w))+1;
	if (k > inExt[2*j+1]) inExt[2*j+1] = k;
	}
      }
    else
      {
      for (j = 0; j < 3; j++) 
	{
	k = int((point[j]+idX*xAxis[j])/w + 0.5);
	if (k < inExt[2*j]) inExt[2*j] = k; 
	k = int((point[j]+idX*xAxis[j])/w + 0.5);
	if (k > inExt[2*j+1]) inExt[2*j+1] = k;
	}
      }
    }
  matrix->Delete();
}


//----------------------------------------------------------------------------
void vtkImageReslice::ExecuteImageInformation() 
{
  int i,j;
  float inPoint[4], outPoint[4];
  float inOrigin[3],maxOut[3],minOut[3],tmp;
  float *inSpacing;
  vtkMatrix4x4 *matrix = vtkMatrix4x4::New();
  int *inExt;

  inExt = this->Input->GetWholeExtent();
  inSpacing = this->Input->GetSpacing();
  this->Input->GetOrigin(inOrigin);
  
  // default extent covers entire input extent
  if (this->OutputExtent[0] == INT_MAX)
    {
    for (i = 0; i < 3; i++)
      {
      minOut[i] = FLT_MAX;
      maxOut[i] = FLT_MIN;
      }
    
    for (i = 0; i < 8; i++)
      {
      inPoint[0] = inOrigin[0] + inExt[i%2]*inSpacing[0];
      inPoint[1] = inOrigin[1] + inExt[2+(i/2)%2]*inSpacing[1];
      inPoint[2] = inOrigin[2] + inExt[4+(i/4)%2]*inSpacing[2];
      inPoint[3] = 1;
      
      if (this->ResliceTransform)
	{
	this->ResliceTransform->GetMatrix(matrix);
	}
      matrix->Invert();
      matrix->MultiplyPoint(inPoint,outPoint);
      
      for (j = 0; j < 3; j++) 
	{
	tmp = outPoint[j]/outPoint[3];
	if (tmp > maxOut[j])
	  maxOut[j] = tmp;
	if (tmp < minOut[j])
	  minOut[j] = tmp;
	}
      }
    
    for (i = 0; i < 3; i++)
      {
      this->OutputExtent[2*i] = 0;
      this->OutputExtent[2*i+1] = 
	int(ceil((maxOut[i]-minOut[i])/this->OutputSpacing[i]));
      }
    
    this->OutputExtent[4] += 1;
    this->OutputExtent[5] += 1;
    }
  // default origin places centre of output over centre of input
  if (this->OutputOrigin[0] == FLT_MAX)
    {
    if (this->ResliceTransform)
      this->ResliceTransform->GetMatrix(matrix);
    matrix->Invert();
    // set translational portion of matrix to zero
    matrix->Element[0][3] = 0;
    matrix->Element[1][3] = 0;
    matrix->Element[2][3] = 0;
    
    for (i = 0; i < 3; i++)
      {
      inPoint[i] = inOrigin[i] + inSpacing[i]*(inExt[2*i]+inExt[2*i+1])*0.5;
      }
    inPoint[3] = 1;
    matrix->MultiplyPoint(inPoint,outPoint);
    for (i = 0; i < 3; i++)
      {
      this->OutputOrigin[i] = outPoint[i]/outPoint[3] 
	- this->OutputSpacing[i]
	*(this->OutputExtent[2*i]+this->OutputExtent[2*i+1])*0.5;
      }
    }
  
  this->Output->SetWholeExtent(this->OutputExtent);
  this->Output->SetSpacing(this->OutputSpacing);
  this->Output->SetOrigin(this->OutputOrigin);
  matrix->Delete();
}

// Do trilinear interpolation of the input data 'inPtr' of extent 'inExt'
// at the 'point'.  The result is placed at 'outPtr'.  
// If the lookup data is beyond the extent 'inExt', set 'outPtr' to
// the background color 'background'.  
// The number of scalar components in the data is 'numscalars'
template <class T>
static void TrilinearWithCheck(float *point, T *inPtr, T *outPtr,
			       T *background, int numscalars, 
			       int inExt[6], int inInc[3])
{
  int i;
  float x,y,z,fx,fy,fz,rx,ry,rz,ryrz,ryfz,fyrz,fyfz;
  int inIdX, inIdY, inIdZ, doInterpX, doInterpY, doInterpZ; 

  inIdX = int(x=point[0])-inExt[0];
  inIdY = int(y=point[1])-inExt[2];
  inIdZ = int(z=point[2])-inExt[4];

  // the doInterpX,Y,Z variables are 0 if interpolation
  // does not have to be done in the specified direction,
  // i.e. if the x, y or z lookup indices have no fractional
  // component. 
  
  doInterpX = (x != int(x));
  doInterpY = (y != int(y));
  doInterpZ = (z != int(z));

  if (inIdX < 0 || inIdX+doInterpX > inExt[1]-inExt[0]
      || inIdY < 0 || inIdY+doInterpY > inExt[3]-inExt[2]
      || inIdZ < 0 || inIdZ+doInterpZ > inExt[5]-inExt[4] )
    {// out of bounds: clear to background color 
    for (i = 0; i < numscalars; i++) 
      *outPtr++ = *background++;
    }
  else 
    {// do trilinear interpolation
    int factX,factY,factZ,factX1,factY1,factZ1;
    T v000,v001,v010,v011,v100,v101,v110,v111;
    
    factX = inIdX*inInc[0];
    factY = inIdY*inInc[1];
    factZ = inIdZ*inInc[2];
    
    factX1 = (inIdX+doInterpX)*inInc[0];
    factY1 = (inIdY+doInterpY)*inInc[1];
    factZ1 = (inIdZ+doInterpZ)*inInc[2];
    
    for (i = 0; i < numscalars; i++) 
      {
      v000 = *(inPtr+factX+factY+factZ);
      v001 = *(inPtr+factX+factY+factZ1);
      v010 = *(inPtr+factX+factY1+factZ);
      v011 = *(inPtr+factX+factY1+factZ1);
      v100 = *(inPtr+factX1+factY+factZ);
      v101 = *(inPtr+factX1+factY+factZ1);
      v110 = *(inPtr+factX1+factY1+factZ);
      v111 = *(inPtr+factX1+factY1+factZ1);
      
      rx = 1 - (fx = x-(inIdX+inExt[0]));
      ry = 1 - (fy = y-(inIdY+inExt[2]));
      rz = 1 - (fz = z-(inIdZ+inExt[4]));
      
      ryrz = ry*rz;
      ryfz = ry*fz;
      fyrz = fy*rz;
      fyfz = fy*fz;
	      
      *outPtr++ = 
	T(rx*(ryrz*v000+ryfz*v001+fyrz*v010+fyfz*v011)
	  + fx*(ryrz*v100+ryfz*v101+fyrz*v110+fyfz*v111)); 
      inPtr++;
      } 
    }
}			  

// Do nearest-neighbor interpolation of the input data 'inPtr' of extent 
// 'inExt' at the 'point'.  The result is placed at 'outPtr'.  
// If the lookup data is beyond the extent 'inExt', set 'outPtr' to
// the background color 'background'.  
// The number of scalar components in the data is 'numscalars'
template <class T>
static void NearestNeighborWithCheck(float *point, T *inPtr, T *outPtr,
				     T *background, int numscalars, 
				     int inExt[6], int inInc[3])
{
  int i;
  int inIdX = int(point[0]+0.5)-inExt[0];
  int inIdY = int(point[1]+0.5)-inExt[2];
  int inIdZ = int(point[2]+0.5)-inExt[4];
  
  if (inIdX < 0 || inIdX > inExt[1]-inExt[0]
      || inIdY < 0 || inIdY > inExt[3]-inExt[2]
      || inIdZ < 0 || inIdZ > inExt[5]-inExt[4] )
    {
    for (i = 0; i < numscalars; i++)
      *outPtr++ = *background++;
    }
  else 
    {
    for (i = 0; i < numscalars; i++)
      *outPtr++ = *(inPtr+inIdX*inInc[0]+inIdY*inInc[1]
		    +inIdZ*inInc[2]);
    }
} 

//----------------------------------------------------------------------------
// This templated function executes the filter for any type of data.
// (this one function is pretty much the be-all and end-all of the
// filter)
template <class T>
static void vtkImageResliceExecute(vtkImageReslice *self,
				   vtkImageData *inData, T *inPtr,
				   vtkImageData *outData, T *outPtr,
				   int outExt[6], int id)
{
  int i, numscalars;
  int idX, idY, idZ;
  int outIncX, outIncY, outIncZ;
  int inIdX, inIdY, inIdZ;
  int inExt[6], inInc[3];
  unsigned long count = 0;
  unsigned long target;
  float inPoint[4],outPoint[4];
  T *background;
  vtkMatrix4x4 *matrix = vtkMatrix4x4::New();
  void (*interpolate)(float *point, T *inPtr, T *outPtr,
		      T *background, int numscalars, 
		      int inExt[6], int inInc[3]);
  
  // find maximum input range
  inData->GetExtent(inExt);
  
  target = (unsigned long)
    ((outExt[5]-outExt[4]+1)*(outExt[3]-outExt[2]+1)/50.0);
  target++;
  
  // Get Increments to march through data 
  inData->GetIncrements(inInc);
  outData->GetContinuousIncrements(outExt, outIncX, outIncY, outIncZ);
  numscalars = inData->GetNumberOfScalarComponents();
  
  // change transform matrix so that instead of taking 
  // input coords -> output coords it takes output indices -> input indices
  self->ComputeIndexMatrix(matrix);
  
  // Color for area outside of input volume extent  
  background = new T[numscalars];
  for (i = 0; i < numscalars; i++)
    background[i] = T(self->GetBackgroundLevel());

  // Set interpolation method
  if (self->GetInterpolate())
    interpolate = &TrilinearWithCheck;
  else
    interpolate = &NearestNeighborWithCheck;

  // Loop through ouput pixels
  for (idZ = outExt[4]; idZ <= outExt[5]; idZ++)
    {
    for (idY = outExt[2]; idY <= outExt[3]; idY++)
      {
      if (!id) 
	{
	if (!(count%target)) 
	  {
	  self->UpdateProgress(count/(50.0*target));
	  }
	count++;
	}
      
      for (idX = outExt[0]; idX <= outExt[1]; idX++)
	{
	outPoint[0] = idX;
	outPoint[1] = idY;
	outPoint[2] = idZ;
	outPoint[3] = 1;

	matrix->MultiplyPoint(outPoint,inPoint); // apply transform

	inPoint[0] /= inPoint[3]; // deal with w if the transform
	inPoint[1] /= inPoint[3]; //   was a perspective transform
	inPoint[2] /= inPoint[3];
	inPoint[3] = 1;
	
	interpolate(inPoint, inPtr, outPtr, background, 
		    numscalars, inExt, inInc);

	outPtr += numscalars; 
	}
      outPtr += outIncY;
      }
    outPtr += outIncZ;
    }
  delete [] background;
  matrix->Delete();
}

//----------------------------------------------------------------------------
// helper functions for vtkOptimizedExecute()

// find approximate intersection of line with the plane x = x_min,
// y = y_min, or z = z_min (lower limit of data extent) 

static int intersectionLow(double *point, float *axis, int *sign,
			   int *limit, int ai)
{
  // approximate value of r
  int r;
  double rd = ((limit[ai]+0.5)*point[3]-point[ai])
    /(axis[ai]-(limit[ai]+0.5)*axis[3]) + 0.5;
  
  if (rd < INT_MIN) 
    r = INT_MIN;
  else if (rd > INT_MAX)
    r = INT_MAX;
  else
    r = int(rd);
  
  // move back and forth to find the point just inside the extent
  while (int( (point[ai]+r*axis[ai])/double(point[3]+r*axis[3]) + 0.5 ) 
	 < limit[ai])
    r += sign[ai];
  
  while (int( (point[ai]+(r-sign[ai])*axis[ai])
	      /double(point[3]+(r-sign[ai])*axis[3]) + 0.5 ) 
	 >= limit[ai])
    r -= sign[ai];
  
  return r;
}

// same as above, but for x = x_max
static int intersectionHigh(double *point, float *axis, int *sign, 
			    int *limit, int ai)
{
  int r;
  double rd = ((limit[ai]-0.5)*point[3]-point[ai])
    /(axis[ai]-(limit[ai]-0.5)*axis[3]) + 0.5; 
    
  if (rd < INT_MIN) 
    r = INT_MIN;
  else if (rd > INT_MAX)
    r = INT_MAX;
  else
    r = int(rd);
  
  while (int( (point[ai]+r*axis[ai])/double(point[3]+r*axis[3]) + 0.5 ) 
	 > limit[ai])
    r -= sign[ai];
  
  while (int( (point[ai]+(r+sign[ai])*axis[ai])
	      /double(point[3]+(r+sign[ai])*axis[3]) + 0.5 ) 
	 <= limit[ai])
    r += sign[ai];
  
  return r;
}

static int isBounded(double *point, float *xAxis, int *inMin, 
		     int *inMax, int ai, int r)
{
  int bi = ai+1; 
  int ci = ai+2;
  if (bi > 2) bi -= 3;  // coordinate index must be 0, 1 or 2
  if (ci > 2) ci -= 3;
  double w = point[3]+r*xAxis[3];
  int bp = int((point[bi]+r*xAxis[bi])/w + 0.5);
  int cp = int((point[ci]+r*xAxis[ci])/w + 0.5);
  
  return (bp >= inMin[bi] && bp <= inMax[bi] &&
	  cp >= inMin[ci] && cp <= inMax[ci]);
}

// this huge mess finds out where the current output raster
// line intersects the input volume 
int vtkImageReslice::FindExtent(int& r1, int& r2, double *point, float *xAxis, 
		      int *inMin, int *inMax)
{
  int i, ix, iy, iz;
  int sign[3];
  int indx1[4],indx2[4];
  double w1,w2;

  // find signs of components of x axis 
  // (this is complicated due to the homogenous coordinate)
  for (i = 0; i < 3; i++)
    {
    if (point[i]/point[3] <= (point[i]+xAxis[i])/(point[3]+xAxis[3]))
      {
      sign[i] = 1;
      }
    else 
      {
      sign[i] = -1;
      }
    } 
  
  // order components of xAxis from largest to smallest
  
  ix = 0;
  for (i = 1; i < 3; i++)
    {
    if (xAxis[i]*xAxis[i] > xAxis[ix]*xAxis[ix])
      {
      ix = i;
      }
    }
  
  iy = ((ix > 1) ? ix-2 : ix+1);
  iz = ((ix > 0) ? ix-1 : ix+2);

  if (xAxis[iz]*xAxis[iz] > xAxis[iy]*xAxis[iy])
    {
    i = iy;
    iy = iz;
    iz = i;
    }

  r1 = intersectionLow(point,xAxis,sign,inMin,ix);
  r2 = intersectionHigh(point,xAxis,sign,inMax,ix);
  
  // find points of intersections
  // first, find w
  w1 = point[3]+r1*xAxis[3];
  w2 = point[3]+r2*xAxis[3];
  
  for (i = 0; i < 3; i++)
    {
      indx1[i] = int((point[i]+r1*xAxis[i])/w1+0.5);
      indx2[i] = int((point[i]+r2*xAxis[i])/w2+0.5);
    }
  if (isBounded(point,xAxis,inMin,inMax,ix,r1))
    { // passed through x face, check opposing face
    if (isBounded(point,xAxis,inMin,inMax,ix,r2))
      {
      return sign[ix];
      }
    
    if (indx2[iy] < inMin[iy])
      { // check y face
      r2 = intersectionLow(point,xAxis,sign,inMin,iy);
      if (isBounded(point,xAxis,inMin,inMax,iy,r2))
	{
	return sign[ix];
	}
      }
    else if (indx2[iy] > inMax[iy])
      { // check other y face
      r2 = intersectionHigh(point,xAxis,sign,inMax,iy);
      if (isBounded(point,xAxis,inMin,inMax,iy,r2))
	{
	return sign[ix];
	}
      }
    
    if (indx2[iz] < inMin[iz])
      { // check z face
      r2 = intersectionLow(point,xAxis,sign,inMin,iz);
      if (isBounded(point,xAxis,inMin,inMax,iz,r2))
	{
	return sign[ix];
	}
      }
    else if (indx2[iz] > inMax[iz])
      { // check other z face
      r2 = intersectionHigh(point,xAxis,sign,inMax,iz);
      if (isBounded(point,xAxis,inMin,inMax,iz,r2))
	{
	return sign[ix];
	}
      }
    }
  
  if (isBounded(point,xAxis,inMin,inMax,ix,r2))
    { // passed through the opposite x face
    if (indx1[iy] < inMin[iy])
	{ // check y face
	r1 = intersectionLow(point,xAxis,sign,inMin,iy);
	if (isBounded(point,xAxis,inMin,inMax,iy,r1))
	  {
	  return sign[ix];
	  }
	}
    else if (indx1[iy] > inMax[iy])
      { // check other y face
      r1 = intersectionHigh(point,xAxis,sign,inMax,iy);
      if (isBounded(point,xAxis,inMin,inMax,iy,r1))
	{
	return sign[ix];
	}
      }
    
    if (indx1[iz] < inMin[iz])
      { // check z face
      r1 = intersectionLow(point,xAxis,sign,inMin,iz);
      if (isBounded(point,xAxis,inMin,inMax,iz,r1))
	{
	return sign[ix];
	}
      }
    else if (indx1[iz] > inMax[iz])
      { // check other z face
      r1 = intersectionHigh(point,xAxis,sign,inMax,iz);
      if (isBounded(point,xAxis,inMin,inMax,iz,r1))
	{
	return sign[ix];
	}
      }
    }
  
  if ((indx1[iy] >= inMin[iy] && indx2[iy] < inMin[iy]) ||
      (indx1[iy] < inMin[iy] && indx2[iy] >= inMin[iy]))
    { // line might pass through bottom face
    r1 = intersectionLow(point,xAxis,sign,inMin,iy);
    if (isBounded(point,xAxis,inMin,inMax,iy,r1))
      {
      if ((indx1[iy] <= inMax[iy] && indx2[iy] > inMax[iy]) ||
	  (indx1[iy] > inMax[iy] && indx2[iy] <= inMax[iy]))
	{ // line might pass through top face
	r2 = intersectionHigh(point,xAxis,sign,inMax,iy);
	if (isBounded(point,xAxis,inMin,inMax,iy,r2))
	  {
	  return sign[iy];
	  }
	}
      
      if (indx1[iz] < inMin[iz] && indx2[iy] < inMin[iy] ||
	  indx2[iz] < inMin[iz] && indx1[iy] < inMin[iy])
	{ // line might pass through in-to-screen face
	r2 = intersectionLow(point,xAxis,sign,inMin,iz);
	if (isBounded(point,xAxis,inMin,inMax,iz,r2))
	  {
	  return sign[iy];
	  }
	}
      else if (indx1[iz] > inMax[iz] && indx2[iy] < inMin[iy] ||
	       indx2[iz] > inMax[iz] && indx1[iy] < inMin[iy])
	{ // line might pass through out-of-screen face
	r2 = intersectionHigh(point,xAxis,sign,inMax,iz);
	if (isBounded(point,xAxis,inMin,inMax,iz,r2))
	  {
	  return sign[iy];
	  }
	} 
      }
    }
  
  if ((indx1[iy] <= inMax[iy] && indx2[iy] > inMax[iy]) ||
      (indx1[iy] > inMax[iy] && indx2[iy] <= inMax[iy]))
    { // line might pass through top face
    r2 = intersectionHigh(point,xAxis,sign,inMax,iy);
    if (isBounded(point,xAxis,inMin,inMax,iy,r2))
      {
      if (indx1[iz] < inMin[iz] && indx2[iy] > inMax[iy] ||
	  indx2[iz] < inMin[iz] && indx1[iy] > inMax[iy])
	{ // line might pass through in-to-screen face
	r1 = intersectionLow(point,xAxis,sign,inMin,iz);
	if (isBounded(point,xAxis,inMin,inMax,iz,r1))
	  {
	  return sign[iy];
	  }
	}
      else if (indx1[iz] > inMax[iz] && indx2[iy] > inMax[iy] || 
	       indx2[iz] > inMax[iz] && indx1[iy] > inMax[iy])
	{ // line might pass through out-of-screen face
	r1 = intersectionHigh(point,xAxis,sign,inMax,iz);
	if (isBounded(point,xAxis,inMin,inMax,iz,r1))
	  {
	  return sign[iy];
	  }
	}
      } 
    }
  
  if ((indx1[iz] >= inMin[iz] && indx2[iz] < inMin[iz]) ||
      (indx1[iz] < inMin[iz] && indx2[iz] >= inMin[iz]))
    { // line might pass through in-to-screen face
    r1 = intersectionLow(point,xAxis,sign,inMin,iz);
    if (isBounded(point,xAxis,inMin,inMax,iz,r1))
      {
      if (indx1[iz] > inMax[iz] || indx2[iz] > inMax[iz])
	{ // line might pass through out-of-screen face
	r2 = intersectionHigh(point,xAxis,sign,inMax,iz);
	if (isBounded(point,xAxis,inMin,inMax,iz,r2))
	  {
	  return sign[iz];
	  }
	}
      }
    }
  
  r1 = r2 = -1;
  return 1;
}


// The vtkOptimizedExecute() function uses an optimization which
// is conceptually simple, but complicated to implement.

// In the un-optimized version, each output voxel
// is converted into a set of look-up indices for the input data;
// then, the indices are checked to ensure they lie within the
// input data extent.

// In the optimized version below, the check is done in reverse:
// it is first determined which output voxels map to look-up indices
// within the input data extent.  Then, further calculations are
// done only for those voxels.  This means that 1) minimal work
// is done for voxels which map to regions outside fo the input
// extent (they are just set to the background color) and 2)
// the inner loops of the look-up and interpolation are
// tightened relative to the un-uptimized version. 

template <class T>
static void vtkOptimizedExecute(vtkImageReslice *self,
				vtkImageData *inData, T *inPtr,
				vtkImageData *outData, T *outPtr,
				int outExt[6], int id)
{
  int i, numscalars;
  int idX, idY, idZ;
  int outIncX, outIncY, outIncZ;
  int inIdX, inIdY, inIdZ;
  int inExt[6];
  int inMax[3], inMin[3];
  int inInc[3];
  unsigned long count = 0;
  unsigned long target;
  int r1,r2,s1,s2;
  double inPoint0[4];
  double inPoint1[4];
  float inPoint[4];
  float *xAxis, *yAxis, *zAxis, *origin;
  T *background;
  vtkMatrix4x4 *matrix = vtkMatrix4x4::New();
  double w;

  // find maximum input range
  inData->GetExtent(inExt);
  
  for (i = 0; i < 3; i++)
    {
    inMin[i] = inExt[2*i];
    inMax[i] = inExt[2*i+1];
    }
  
  target = (unsigned long)
    ((outExt[5]-outExt[4]+1)*(outExt[3]-outExt[2]+1)/50.0);
  target++;
  
  // Get Increments to march through data 
  inData->GetIncrements(inInc);
  outData->GetContinuousIncrements(outExt, outIncX, outIncY, outIncZ);
  numscalars = inData->GetNumberOfScalarComponents();

  // set up background levels
  background = new T[numscalars];
  for (i = 0; i < numscalars; i++)
    background[i] = T(self->GetBackgroundLevel());

  // change transform matrix so that instead of taking 
  // input coords -> output coords it takes output indices -> input indices
  self->ComputeIndexMatrix(matrix);
  
  // break matrix into a set of axes plus an origin
  // (this allows us to calculate the transform Incrementally)
  matrix->Transpose();
  xAxis = matrix->Element[0];
  yAxis = matrix->Element[1];
  zAxis = matrix->Element[2];
  origin = matrix->Element[3];
  
  // Loop through output pixels
  for (idZ = outExt[4]; idZ <= outExt[5]; idZ++)
    {
    inPoint0[0] = origin[0]+idZ*zAxis[0]; // incremental transform
    inPoint0[1] = origin[1]+idZ*zAxis[1]; 
    inPoint0[2] = origin[2]+idZ*zAxis[2]; 
    inPoint0[3] = origin[3]+idZ*zAxis[3]; 
    
    for (idY = outExt[2]; idY <= outExt[3]; idY++)
      {
      inPoint1[0] = inPoint0[0]+idY*yAxis[0]; // incremental transform
      inPoint1[1] = inPoint0[1]+idY*yAxis[1];
      inPoint1[2] = inPoint0[2]+idY*yAxis[2];
      inPoint1[3] = inPoint0[3]+idY*yAxis[3];
      
      if (!id) 
	{
	if (!(count%target)) 
	  {
	  self->UpdateProgress(count/(50.0*target));
	  }
	count++;
	}
      
      // find intersections of x raster line with the input extent
      if (self->FindExtent(r1,r2,inPoint1,xAxis,inMin,inMax) < 0)
	{
	i = r1;
	r1 = r2;
	r2 = i;
	}

      // bound r1,r2 within reasonable limits
      if (r1 < outExt[0]) r1 = outExt[0];
      if (r2 > outExt[1]) r2 = outExt[1];
      if (r1 > r2) r1 = r2+1;
      if (r2 < r1) r2 = r1-1;

      // clear pixels to left of input extent
      if (numscalars == 1) // optimize for single scalar
	for (idX = outExt[0]; idX < r1; idX++) 
	  *outPtr++ = background[0];
      else             // multiple scalars
	for (idX = outExt[0]; idX < r1; idX++)
	  for (i = 0; i < numscalars; i++)
	    *outPtr++ = background[i];
      
      // interpolate pixels while inside input extent
      if (self->GetInterpolate())
	{ // do trilinear interpolation
	float x,y,z,fx,fy,fz,rx,ry,rz,ryrz,ryfz,fyrz,fyfz;
	T v000,v001,v010,v011,v100,v101,v110,v111;
	int factX,factY,factZ,factX1,factY1,factZ1;
	int doInterpX, doInterpY, doInterpZ;

	// check out the questionable region

	// adjust left bound if necessary
	for (s1 = r1; s1 <= r2; s1++)
	  {	  
	  w = inPoint1[3]+idX*xAxis[3]; // don't forget w!  
	  inIdX = int(x=float((inPoint1[0]+s1*xAxis[0])/w))-inExt[0];
	  inIdY = int(y=float((inPoint1[1]+s1*xAxis[1])/w))-inExt[2];
	  inIdZ = int(z=float((inPoint1[2]+s1*xAxis[2])/w))-inExt[4];

	  if (inIdX >= 0 && inIdX+1 <= inExt[1]-inExt[0]
	      && inIdY >= 0 && inIdY+1 <= inExt[3]-inExt[2]
	      && inIdZ >= 0 && inIdZ+1 <= inExt[5]-inExt[4] )
	    break;
	  }
	// adjust right bound if necessary
	for (s2 = r2; s2 >= s1; s2--)
	  {	  
	  w = inPoint1[3]+idX*xAxis[3]; // don't forget w!  
	  inIdX = int(x=float((inPoint1[0]+s2*xAxis[0])/w))-inExt[0];
	  inIdY = int(y=float((inPoint1[1]+s2*xAxis[1])/w))-inExt[2];
	  inIdZ = int(z=float((inPoint1[2]+s2*xAxis[2])/w))-inExt[4];

	  if (inIdX >= 0 && inIdX+1 <= inExt[1]-inExt[0]
	      && inIdY >= 0 && inIdY+1 <= inExt[3]-inExt[2]
	      && inIdZ >= 0 && inIdZ+1 <= inExt[5]-inExt[4] )
	    break;
	  }

	// go through the left part of the questionable region

	for (; idX < s1; idX++)
	  {
	  w = inPoint1[3]+idX*xAxis[3];
	  inPoint[0] = (inPoint1[0]+idX*xAxis[0])/w;
	  inPoint[1] = (inPoint1[1]+idX*xAxis[1])/w;
	  inPoint[2] = (inPoint1[2]+idX*xAxis[2])/w;
	  inPoint[3] = 1;

	  TrilinearWithCheck(inPoint, inPtr, outPtr, background, 
			     numscalars, inExt, inInc);

	  outPtr += numscalars;
	  }
	
	// go through the contiguous section, no checks are
	// required
  
	for (; idX <= s2; idX++)
	  {
	  w = inPoint1[3]+idX*xAxis[3]; // don't forget w!  
	  inIdX = int(x=float((inPoint1[0]+idX*xAxis[0])/w))-inExt[0];
	  inIdY = int(y=float((inPoint1[1]+idX*xAxis[1])/w))-inExt[2];
	  inIdZ = int(z=float((inPoint1[2]+idX*xAxis[2])/w))-inExt[4];
	  factX = inIdX*inInc[0];
	  factY = inIdY*inInc[1];
	  factZ = inIdZ*inInc[2];
	  factX1 = (inIdX+1)*inInc[0];
	  factY1 = (inIdY+1)*inInc[1];
	  factZ1 = (inIdZ+1)*inInc[2];
	  
	  for (i = 0; i < numscalars; i++) 
	    {
	    v000 = *(inPtr+factX+factY+factZ);
	    v001 = *(inPtr+factX+factY+factZ1);
	    v010 = *(inPtr+factX+factY1+factZ);
	    v011 = *(inPtr+factX+factY1+factZ1);
	    v100 = *(inPtr+factX1+factY+factZ);
	    v101 = *(inPtr+factX1+factY+factZ1);
	    v110 = *(inPtr+factX1+factY1+factZ);
	    v111 = *(inPtr+factX1+factY1+factZ1);
	    
	    rx = 1 - (fx = x-(inIdX+inExt[0]));
	    ry = 1 - (fy = y-(inIdY+inExt[2]));
	    rz = 1 - (fz = z-(inIdZ+inExt[4]));
	    
	    ryrz = ry*rz;
	    ryfz = ry*fz;
	    fyrz = fy*rz;
	    fyfz = fy*fz;
	    
	    *outPtr++ = 
	      T(rx*(ryrz*v000+ryfz*v001+fyrz*v010+fyfz*v011)
		+ fx*(ryrz*v100+ryfz*v101+fyrz*v110+fyfz*v111));
	    
	    inPtr++;
	    }
	  inPtr -= numscalars;
	  }

	// go through the left part of the questionable region

	for (; idX <= r2; idX++)
	  {
	  w = inPoint1[3]+idX*xAxis[3];
	  inPoint[0] = (inPoint1[0]+idX*xAxis[0])/w;
	  inPoint[1] = (inPoint1[1]+idX*xAxis[1])/w;
	  inPoint[2] = (inPoint1[2]+idX*xAxis[2])/w;
	  inPoint[3] = 1;

	  TrilinearWithCheck(inPoint, inPtr, outPtr, background, 
			     numscalars, inExt, inInc);

	  outPtr += numscalars;
	  }
	
	}
      else
	{ // do nearest-neighbor interpolation
	if (numscalars == 1)
	  { // optimize for single scalar (keep that inner loop tight!)
	  for (; idX <= r2; idX++)
	    {
	    w = inPoint1[3]+idX*xAxis[3]; // don't forget w!  
	    inIdX = int((inPoint1[0]+idX*xAxis[0])/w+0.5)-inExt[0];
	    inIdY = int((inPoint1[1]+idX*xAxis[1])/w+0.5)-inExt[2];
	    inIdZ = int((inPoint1[2]+idX*xAxis[2])/w+0.5)-inExt[4];
	    
	    *outPtr++ = *(inPtr+inIdX*inInc[0]+inIdY*inInc[1]
			  +inIdZ*inInc[2]);
	    }
	  }
	else
	  { // multiple scalar values
	  for (; idX <= r2; idX++)
	    {
	    w = inPoint1[3]+idX*xAxis[3]; // don't forget w!  
	    inIdX = int((inPoint1[0]+idX*xAxis[0])/w+0.5)-inExt[0];
	    inIdY = int((inPoint1[1]+idX*xAxis[1])/w+0.5)-inExt[2];
	    inIdZ = int((inPoint1[2]+idX*xAxis[2])/w+0.5)-inExt[4];
	    
	    for (i = 0; i < numscalars; i++)
	      *outPtr++ = *(inPtr+inIdX*inInc[0]+inIdY*inInc[1]
			    +inIdZ*inInc[2]);
	    }
	  }
	}
      
      // clear pixels to right of input extent
      if (numscalars == 1) // optimize for single scalar
	for (; idX <= outExt[1]; idX++) 
	  *outPtr++ = background[0];
      else // multiple scalars
	for (; idX <= outExt[1]; idX++)
	  for (i = 0; i < numscalars; i++)
	    *outPtr++ = background[i];
      
      outPtr += outIncY;
      }
    outPtr += outIncZ;
    }
  matrix->Delete();
  delete [] background;
}

//----------------------------------------------------------------------------
// This method is passed a input and output region, and executes the filter
// algorithm to fill the output from the input.
// It just executes a switch statement to call the correct function for
// the regions data types.
void vtkImageReslice::ThreadedExecute(vtkImageData *inData, 
				      vtkImageData *outData,
				      int outExt[6], int id)
{
  void *inPtr = inData->GetScalarPointerForExtent(inData->GetExtent());
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
  
  if (this->Optimization)
    {
    switch (inData->GetScalarType())
      {
      case VTK_FLOAT:
	vtkOptimizedExecute(this, inData, (float *)(inPtr), 
			       outData, (float *)(outPtr),outExt, id);
	break;
      case VTK_INT:
	vtkOptimizedExecute(this, inData, (int *)(inPtr), 
			       outData, (int *)(outPtr),outExt, id);
	break;
      case VTK_SHORT:
	vtkOptimizedExecute(this, inData, (short *)(inPtr), 
			       outData, (short *)(outPtr),outExt, id);
	break;
      case VTK_UNSIGNED_SHORT:
	vtkOptimizedExecute(this, inData, (unsigned short *)(inPtr), 
			       outData, (unsigned short *)(outPtr),outExt,id);
	break;
      case VTK_UNSIGNED_CHAR:
	vtkOptimizedExecute(this, inData, (unsigned char *)(inPtr), 
			       outData, (unsigned char *)(outPtr),outExt, id);
	break;
      default:
	vtkErrorMacro(<< "Execute: Unknown input ScalarType");
	return;
      }
    }
  else
    {
    switch (inData->GetScalarType())
      {
      case VTK_FLOAT:
	vtkImageResliceExecute(this, inData, (float *)(inPtr), 
			       outData, (float *)(outPtr),outExt, id);
	break;
      case VTK_INT:
	vtkImageResliceExecute(this, inData, (int *)(inPtr), 
			       outData, (int *)(outPtr),outExt, id);
	break;
      case VTK_SHORT:
	vtkImageResliceExecute(this, inData, (short *)(inPtr), 
			       outData, (short *)(outPtr),outExt, id);
	break;
      case VTK_UNSIGNED_SHORT:
	vtkImageResliceExecute(this, inData, (unsigned short *)(inPtr), 
			       outData, (unsigned short *)(outPtr),outExt,id);
	break;
      case VTK_UNSIGNED_CHAR:
	vtkImageResliceExecute(this, inData, (unsigned char *)(inPtr), 
			       outData, (unsigned char *)(outPtr),outExt, id);
	break;
      default:
	vtkErrorMacro(<< "Execute: Unknown input ScalarType");
	return;
      }
    }
}


