/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkImageReslice.cxx
  Language:  C++
  Date:      $Date$

  Version:   $Revision$
  Thanks:    Thanks to David G Gobbi who developed this class.

Copyright (c) 1993-2001 Ken Martin, Will Schroeder, Bill Lorensen 
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.

 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

 * Neither name of Ken Martin, Will Schroeder, or Bill Lorensen nor the names
   of any contributors may be used to endorse or promote products derived
   from this software without specific prior written permission.

 * Modified source versions must be plainly marked as such, and must not be
   misrepresented as being the original software.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS ``AS IS''
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=========================================================================*/
#include <limits.h>
#include <float.h>
#include <math.h>
#include "vtkImageReslice.h"
#include "vtkMath.h"
#include "vtkTransform.h"
#include "vtkObjectFactory.h"

//----------------------------------------------------------------------------
vtkImageReslice* vtkImageReslice::New()
{
  // First try to create the object from the vtkObjectFactory
  vtkObject* ret = vtkObjectFactory::CreateInstance("vtkImageReslice");
  if(ret)
    {
    return (vtkImageReslice*)ret;
    }
  // If the factory was unable to create the object, then create it here.
  return new vtkImageReslice;
}

//----------------------------------------------------------------------------
vtkImageReslice::vtkImageReslice()
{
  this->OutputSpacing[0] = 1;
  this->OutputSpacing[1] = 1;
  this->OutputSpacing[2] = 1;

  // flag to set defaults later
  this->OutputOrigin[0] = FLT_MAX;
  this->OutputOrigin[1] = FLT_MAX;
  this->OutputOrigin[2] = FLT_MAX;

  // ditto
  this->OutputExtent[0] = this->OutputExtent[1] = INT_MAX;
  this->OutputExtent[2] = this->OutputExtent[3] = INT_MAX;
  this->OutputExtent[4] = this->OutputExtent[5] = INT_MAX;

  this->OutputAlwaysCenteredOnInput = 0;
  
  this->Wrap = 0; // don't wrap
  this->Mirror = 0; // don't mirror
  this->InterpolationMode = VTK_RESLICE_NEAREST; // no interpolation
  this->Optimization = 1; // optimizations seem to finally be stable...

  this->BackgroundColor[0] = 0;
  this->BackgroundColor[1] = 0;
  this->BackgroundColor[2] = 0;
  this->BackgroundColor[3] = 0;

  this->ResliceAxes = NULL;
  this->ResliceTransform = NULL;
  this->IndexMatrix = NULL;
}


//----------------------------------------------------------------------------

vtkImageReslice::~vtkImageReslice()
{
  this->SetResliceTransform(NULL);
  this->SetResliceAxes(NULL);
  if (this->IndexMatrix)
    {
    this->IndexMatrix->Delete();
    }
}

//----------------------------------------------------------------------------
void vtkImageReslice::PrintSelf(ostream& os, vtkIndent indent)
{
  vtkImageToImageFilter::PrintSelf(os,indent);

  os << indent << "ResliceAxes: " << this->ResliceAxes << "\n";
  if (this->ResliceAxes)
    {
    this->ResliceAxes->PrintSelf(os,indent.GetNextIndent());
    }
  os << indent << "ResliceTransform: " << this->ResliceTransform << "\n";
  if (this->ResliceTransform)
    {
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
  os << indent << "OutputAlwaysCenteredOnInput: " << 
    (this->OutputAlwaysCenteredOnInput ? "On\n":"Off\n");
  os << indent << "Wrap: " << (this->Wrap ? "On\n":"Off\n");
  os << indent << "Mirror: " << (this->Mirror ? "On\n":"Off\n");
  os << indent << "InterpolationMode: " 
     << this->GetInterpolationModeAsString() << "\n";
  os << indent << "Optimization: " << (this->Optimization ? "On\n":"Off\n");
  os << indent << "BackgroundColor: " << this->BackgroundColor[0] << " " <<
    this->BackgroundColor[1] << " " << this->BackgroundColor[2] << " " <<
    this->BackgroundColor[3] << "\n";
}

//----------------------------------------------------------------------------
// Account for the MTime of the transform and its matrix when determinging
// the MTime of the filter

unsigned long int vtkImageReslice::GetMTime()
{
  unsigned long mTime=this->vtkObject::GetMTime();
  unsigned long time;

  if ( this->ResliceTransform != NULL )
    {
    time = this->ResliceTransform->GetMTime();
    mTime = ( time > mTime ? time : mTime );
    if (this->ResliceTransform->IsA("vtkHomogeneousTransform"))
      {
      time = ((vtkHomogeneousTransform *)this->ResliceTransform)
	->GetMatrix()->GetMTime();
      mTime = ( time > mTime ? time : mTime );
      }    
    }
  if ( this->ResliceAxes != NULL)
    {
    time = this->ResliceAxes->GetMTime();
    mTime = ( time > mTime ? time : mTime );
    }

  return mTime;
}

//----------------------------------------------------------------------------
// fast floor() function for converting a float to an int
// (the floor() implementation on some computers is much slower than this,
// because they require some 'exact' behaviour that we don't).

static inline int vtkResliceFloor(float x, float &f)
{
  int ix = int(x);
  f = x-ix;
  if (f < 0) { f = x - (--ix); }

  return ix;
}

static inline int vtkResliceFloor(float x)
{
  int ix = int(x);
  if (x-ix < 0) { ix--; }

  return ix;
}

static inline int vtkResliceCeil(float x)
{
  int ix = int(x);
  if (x-ix > 0) { ix++; }

  return ix;
}

//----------------------------------------------------------------------------
void vtkImageReslice::ComputeInputUpdateExtent(int inExt[6], 
					       int outExt[6])
{
  if (this->ResliceTransform)
    {
    this->ResliceTransform->Update();
    if (!this->ResliceTransform->IsA("vtkHomogeneousTransform"))
      { // set the input to the whole extent if the transform
	// is nonlinear
      this->GetInput()->GetWholeExtent(inExt);
      return;
      }
    }

  if (this->Optimization)
    {
    this->OptimizedComputeInputUpdateExtent(inExt,outExt);
    return;
    }

  int i,j,k;
  float point[4],f;
  float *inSpacing,*inOrigin,*outSpacing,*outOrigin,inInvSpacing[3];

  int wrap = (this->GetWrap() || this->GetInterpolationMode() != VTK_RESLICE_NEAREST);
  
  inOrigin = this->GetInput()->GetOrigin();
  inSpacing = this->GetInput()->GetSpacing();
  outOrigin = this->GetOutputOrigin();
  outSpacing = this->GetOutputSpacing();

  // save effor later: invert inSpacing
  inInvSpacing[0] = 1.0f/inSpacing[0];
  inInvSpacing[1] = 1.0f/inSpacing[1];
  inInvSpacing[2] = 1.0f/inSpacing[2];

  for (i = 0; i < 3; i++)
    {
    inExt[2*i] = INT_MAX;
    inExt[2*i+1] = INT_MIN;
    }

  // check the coordinates of the 8 corners of the output extent
  for (i = 0; i < 8; i++)  
    {
    // get output coords
    point[0] = outExt[i%2];
    point[1] = outExt[2+(i/2)%2];
    point[2] = outExt[4+(i/4)%2];

    point[0] = point[0]*outSpacing[0] + outOrigin[0];
    point[1] = point[1]*outSpacing[1] + outOrigin[1];
    point[2] = point[2]*outSpacing[2] + outOrigin[2];
    
    if (this->ResliceAxes)
      {
      point[3] = 1.0f;
      this->ResliceAxes->MultiplyPoint(point,point);
      f = 1.0f/point[3];
      point[0] *= f;
      point[1] *= f;
      point[2] *= f;
      }
    if (this->ResliceTransform)
      {
      this->ResliceTransform->TransformPoint(point,point);
      }  

    point[0] = (point[0] - inOrigin[0])*inInvSpacing[0];
    point[1] = (point[1] - inOrigin[1])*inInvSpacing[1];
    point[2] = (point[2] - inOrigin[2])*inInvSpacing[2];

    // set the extent appropriately according to the interpolation mode 
    if (this->GetInterpolationMode() != VTK_RESLICE_NEAREST)
      {
      int extra = (this->GetInterpolationMode() == VTK_RESLICE_CUBIC); 
      for (j = 0; j < 3; j++) 
	{
	k = vtkResliceFloor(point[j])-extra;
	if (k < inExt[2*j]) 
	  {
	  inExt[2*j] = k;
	  }
	if (wrap)
	  {
	  k = vtkResliceFloor(point[j])+1+extra;	
	  }
	else
	  {
	  k = vtkResliceCeil(point[j])+extra;
	  }	
	if (k > inExt[2*j+1])
	  {
	  inExt[2*j+1] = k;
	  }
	}
      }
    else
      {
      for (j = 0; j < 3; j++) 
	{
	k = vtkResliceFloor(point[j] + 0.5f);
	if (k < inExt[2*j])
	  { 
	  inExt[2*j] = k;
	  } 
	if (k > inExt[2*j+1]) 
	  {
	  inExt[2*j+1] = k;
	  }
	}
      }
    }

  int *wholeExtent = this->GetInput()->GetWholeExtent();
  // Clip, just to make sure we hit _some_ of the input extent
  for (i = 0; i < 3; i++)
    {
    if (inExt[2*i] < wholeExtent[2*i])
      {
      inExt[2*i] = wholeExtent[2*i];
      if (wrap)
	{
	inExt[2*i+1] = wholeExtent[2*i+1];
	}
      }
    if (inExt[2*i+1] > wholeExtent[2*i+1])
      {
      inExt[2*i+1] = wholeExtent[2*i+1];
      if (wrap)
	{
	inExt[2*i] = wholeExtent[2*i];
	}
      }
    if (inExt[2*i] > wholeExtent[2*i+1])
      {
      inExt[2*i] = wholeExtent[2*i+1];
      }
    if (inExt[2*i+1] < wholeExtent[2*i])
      {
      inExt[2*i+1] = wholeExtent[2*i];
      }
    }
}

//----------------------------------------------------------------------------
void vtkImageReslice::ExecuteInformation(vtkImageData *input, 
					 vtkImageData *output) 
{
  int i,j;
  float inPoint[4], outPoint[4];
  float inOrigin[3],maxOut[3],minOut[3],f;
  float *inSpacing;

  int *inWholeExt;

  input->UpdateInformation();
  inWholeExt = input->GetWholeExtent();
  inSpacing = input->GetSpacing();
  input->GetOrigin(inOrigin);
  
  vtkMatrix4x4 *matrix = vtkMatrix4x4::New();

  if (this->ResliceAxes)
    {
    matrix->DeepCopy(this->ResliceAxes);
    }
  if (this->ResliceTransform && 
      this->ResliceTransform->IsA("vtkHomogeneousTransform"))
    {
    vtkMatrix4x4 *transformMatrix = 
      ((vtkHomogeneousTransform *)this->ResliceTransform)->GetMatrix();
    this->ResliceTransform->Update();
    vtkMatrix4x4::Multiply4x4(transformMatrix,matrix,matrix);
    }
  
  // because vtkMatrix4x4::Inverse() doesn't cut it,
  // use vtkMath::InvertMatrix()
  double mat1data[4][4];
  double mat2data[4][4];
  double *mat1[4];
  double *mat2[4];
    
  int tmpIntSpace[4];
  double tmpDoubleSpace[4];

  for (i = 0; i < 4; i++)
    {
    mat1[i] = mat1data[i];
    mat2[i] = mat2data[i];
    for (j = 0; j < 4; j++)
      { 
      mat1[i][j] = matrix->GetElement(i,j);
      }
    }
 
  if (vtkMath::InvertMatrix(mat1,mat2,4,tmpIntSpace,tmpDoubleSpace) == 0)
    {
    vtkErrorMacro(<< "ExecuteInformation: reslicing transform not \
invertible");
    }

  for (i = 0; i < 4; i++)
    {
    for (j = 0; j < 4; j++)
      {
      matrix->SetElement(i,j,mat2[i][j]);
      }
    }

  // default extent covers entire input extent
  if ( this->OutputAlwaysCenteredOnInput || this->OutputExtent[0] == INT_MAX)
    {
    for (i = 0; i < 3; i++)
      {
      minOut[i] = FLT_MAX;
      maxOut[i] = -FLT_MAX;
      }
    
    for (i = 0; i < 8; i++)
      {
      inPoint[0] = inOrigin[0] + inWholeExt[i%2]*inSpacing[0];
      inPoint[1] = inOrigin[1] + inWholeExt[2+(i/2)%2]*inSpacing[1];
      inPoint[2] = inOrigin[2] + inWholeExt[4+(i/4)%2]*inSpacing[2];
      inPoint[3] = 1.0f;
      
      matrix->MultiplyPoint(inPoint,outPoint);
      
      f = 1.0f/outPoint[3];
      outPoint[0] *= f; 
      outPoint[1] *= f; 
      outPoint[2] *= f;

      for (j = 0; j < 3; j++) 
	{
	if (outPoint[j] > maxOut[j])
	  {
	  maxOut[j] = outPoint[j];
	  }
	if (outPoint[j] < minOut[j])
	  {
	  minOut[j] = outPoint[j];
	  }
	}
      }
    
    for (i = 0; i < 3; i++)
      {
      float spacing = this->OutputSpacing[i];
      if (spacing < 0)
	{
	float tmp = maxOut[i];
	maxOut[i] = minOut[i];
	minOut[i] = tmp;
	}
      this->OutputExtent[2*i] = inWholeExt[2*i];
      this->OutputExtent[2*i+1] = inWholeExt[2*i]+
	vtkResliceCeil((maxOut[i]-minOut[i])/spacing);

      if (this->OutputAlwaysCenteredOnInput || this->OutputOrigin[i] == FLT_MAX )
	{
	this->OutputOrigin[i] = minOut[i]-this->OutputExtent[2*i]*spacing;
	}
      }    
    }

  // default origin places centre of output over centre of input
  if (this->OutputOrigin[0] == FLT_MAX )
    {
    for (i = 0; i < 3; i++)
      {
      inPoint[i] = inOrigin[i] + 
	inSpacing[i]*(inWholeExt[2*i]+inWholeExt[2*i+1])*0.5f;
      }
    inPoint[3] = 1.0f;

    matrix->MultiplyPoint(inPoint,outPoint);

    f = 1.0f/outPoint[3];
    for (i = 0; i < 3; i++)
      {
      outPoint[i] *= f;
      this->OutputOrigin[i] = outPoint[i] - this->OutputSpacing[i]
	*(this->OutputExtent[2*i]+this->OutputExtent[2*i+1])*0.5f;
      }
    }
  
  output->SetWholeExtent(this->OutputExtent);
  output->SetSpacing(this->OutputSpacing);
  output->SetOrigin(this->OutputOrigin);
  output->SetScalarType(input->GetScalarType());
  output->SetNumberOfScalarComponents(input->GetNumberOfScalarComponents());

  matrix->Delete();
}

//----------------------------------------------------------------------------
//  Interpolation subroutines and associated code
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
// rounding functions, split and optimized for each type
// (because we don't want to round if the result is a float!)

// in the case of a tie between integers, the larger integer wins.

static inline void vtkResliceRound(float val, unsigned char& rnd)
{
  rnd = (unsigned char)(val+0.5f);
}

static inline void vtkResliceRound(float val, short& rnd)
{
  rnd = (short)((int)(val+32768.5f)-32768);
}

static inline void vtkResliceRound(float val, unsigned short& rnd)
{
  rnd = (unsigned short)(val+0.5f);
}

static inline void vtkResliceRound(float val, int& rnd)
{
  rnd = (int)(floor(val+0.5f));
}

static inline void vtkResliceRound(float val, float& rnd)
{
  rnd = (float)(val);
}

//----------------------------------------------------------------------------
// clamping functions for each type

static inline void vtkResliceClamp(float val, unsigned char& clamp)
{
  if (val < VTK_UNSIGNED_CHAR_MIN)
    { 
    val = VTK_UNSIGNED_CHAR_MIN;
    }
  if (val > VTK_UNSIGNED_CHAR_MAX)
    { 
    val = VTK_UNSIGNED_CHAR_MAX;
    }
  vtkResliceRound(val,clamp);
}

static inline void vtkResliceClamp(float val, short& clamp)
{
  if (val < VTK_SHORT_MIN)
    { 
    val = VTK_SHORT_MIN;
    }
  if (val > VTK_SHORT_MAX)
    { 
    val = VTK_SHORT_MAX;
    }
  vtkResliceRound(val,clamp);
}

static inline void vtkResliceClamp(float val, unsigned short& clamp)
{
  if (val < VTK_UNSIGNED_SHORT_MIN)
    { 
    val = VTK_UNSIGNED_SHORT_MIN;
    }
  if (val > VTK_UNSIGNED_SHORT_MAX)
    { 
    val = VTK_UNSIGNED_SHORT_MAX;
    }
  vtkResliceRound(val,clamp);
}

static inline void vtkResliceClamp(float val, int& clamp)
{
  if (val < VTK_INT_MIN) 
    {
    val = VTK_INT_MIN;
    }
  if (val > VTK_INT_MAX) 
    {
    val = VTK_INT_MAX;
    }
  vtkResliceRound(val,clamp);
}

static inline void vtkResliceClamp(float val, float& clamp)
{
  if (val < VTK_FLOAT_MIN)
    { 
    val = VTK_FLOAT_MIN;
    }
  if (val > VTK_FLOAT_MAX) 
    {
    val = VTK_FLOAT_MAX;
    }
  vtkResliceRound(val,clamp);
}

//----------------------------------------------------------------------------
// copy a pixel, advance the output pointer but not the input pointer

template<class T>
static inline void vtkCopyPixel(T *&out, T *in, int numscalars)
{
  do
    {
    *out++ = *in++;
    }
  while (--numscalars);
}

//----------------------------------------------------------------------------
// Perform a wrap to limit an index to [0,range).
// Ensures correct behaviour when the index is negative.
 
static inline int vtkInterpolateWrap(int num, int range)
{
  if ((num %= range) < 0)
    {
    num += range; // required for some % implementations
    } 
  return num;
}

//----------------------------------------------------------------------------
// Perform a mirror to limit an index to [0,range).
 
static inline int vtkInterpolateMirror(int num, int range)
{
  if (num < 0)
    {
    num = -num-1;
    }
  int count = num/range;
  num %= range;
  if (count & 0x1)
    {
    num = range-num-1;
    }
  return num;
}

//----------------------------------------------------------------------------
// Do trilinear interpolation of the input data 'inPtr' of extent 'inExt'
// at the 'point'.  The result is placed at 'outPtr'.  
// If the lookup data is beyond the extent 'inExt', set 'outPtr' to
// the background color 'background'.  
// The number of scalar components in the data is 'numscalars'
template <class T>
static int vtkTrilinearInterpolation(float *point, T *inPtr, T *outPtr,
				     T *background, int numscalars, 
				     int inExt[6], int inInc[3])
{
  float fx,fy,fz;
  int floorX = vtkResliceFloor(point[0],fx);
  int floorY = vtkResliceFloor(point[1],fy);
  int floorZ = vtkResliceFloor(point[2],fz);

  int inIdX0 = floorX-inExt[0];
  int inIdY0 = floorY-inExt[2];
  int inIdZ0 = floorZ-inExt[4];

  int inIdX1 = inIdX0 + (fx != 0);
  int inIdY1 = inIdY0 + (fy != 0);
  int inIdZ1 = inIdZ0 + (fz != 0);
  
  if (inIdX0 < 0 || inIdX1 > inExt[1]-inExt[0]
      || inIdY0 < 0 || inIdY1 > inExt[3]-inExt[2]
      || inIdZ0 < 0 || inIdZ1 > inExt[5]-inExt[4] )
    {// out of bounds: clear to background color 
    if (background)
      {
      vtkCopyPixel(outPtr,background,numscalars);
      }
    return 0;
    }
  else 
    {// do trilinear interpolation
    int factX = inIdX0*inInc[0];
    int factY = inIdY0*inInc[1];
    int factZ = inIdZ0*inInc[2];

    int factX1 = inIdX1*inInc[0];
    int factY1 = inIdY1*inInc[1];
    int factZ1 = inIdZ1*inInc[2];
    
    int i000 = factX+factY+factZ;
    int i001 = factX+factY+factZ1;
    int i010 = factX+factY1+factZ;
    int i011 = factX+factY1+factZ1;
    int i100 = factX1+factY+factZ;
    int i101 = factX1+factY+factZ1;
    int i110 = factX1+factY1+factZ;
    int i111 = factX1+factY1+factZ1;

    float rx = 1.0f - fx;
    float ry = 1.0f - fy;
    float rz = 1.0f - fz;
      
    float ryrz = ry*rz;
    float ryfz = ry*fz;
    float fyrz = fy*rz;
    float fyfz = fy*fz;

    do
      {
      vtkResliceRound((rx*(ryrz*inPtr[i000]+ryfz*inPtr[i001]+
			   fyrz*inPtr[i010]+fyfz*inPtr[i011])
		       + fx*(ryrz*inPtr[i100]+ryfz*inPtr[i101]+
			     fyrz*inPtr[i110]+fyfz*inPtr[i111])),
		      *outPtr++);
      inPtr++;
      }
    while (--numscalars);

    return 1;
    }
}			  

// trilinear interpolation with wrap-around behaviour
template <class T>
static int vtkTrilinearInterpolationRepeat(float *point, T *inPtr, T *outPtr,
					   T *mirror, int numscalars, 
					   int inExt[6], int inInc[3])
{
  float fx,fy,fz;
  int floorX = vtkResliceFloor(point[0],fx);
  int floorY = vtkResliceFloor(point[1],fy);
  int floorZ = vtkResliceFloor(point[2],fz);

  int inIdX = floorX-inExt[0];
  int inIdY = floorY-inExt[2];
  int inIdZ = floorZ-inExt[4];

  int inExtX = inExt[1]-inExt[0]+1;
  int inExtY = inExt[3]-inExt[2]+1;
  int inExtZ = inExt[5]-inExt[4]+1;

  int factX, factY, factZ;
  int factX1, factY1, factZ1;

  if (mirror)
    {
    factX = vtkInterpolateMirror(inIdX,inExtX)*inInc[0];
    factY = vtkInterpolateMirror(inIdY,inExtY)*inInc[1];
    factZ = vtkInterpolateMirror(inIdZ,inExtZ)*inInc[2];

    factX1 = vtkInterpolateMirror(inIdX+1,inExtX)*inInc[0];
    factY1 = vtkInterpolateMirror(inIdY+1,inExtY)*inInc[1];
    factZ1 = vtkInterpolateMirror(inIdZ+1,inExtZ)*inInc[2];
    }
  else
    {
    factX = vtkInterpolateWrap(inIdX,inExtX)*inInc[0];
    factY = vtkInterpolateWrap(inIdY,inExtY)*inInc[1];
    factZ = vtkInterpolateWrap(inIdZ,inExtZ)*inInc[2];

    factX1 = vtkInterpolateWrap(inIdX+1,inExtX)*inInc[0];
    factY1 = vtkInterpolateWrap(inIdY+1,inExtY)*inInc[1];
    factZ1 = vtkInterpolateWrap(inIdZ+1,inExtZ)*inInc[2];
    }

  int i000 = factX+factY+factZ;
  int i001 = factX+factY+factZ1;
  int i010 = factX+factY1+factZ;
  int i011 = factX+factY1+factZ1;
  int i100 = factX1+factY+factZ;
  int i101 = factX1+factY+factZ1;
  int i110 = factX1+factY1+factZ;
  int i111 = factX1+factY1+factZ1;

  float rx = 1.0f - fx;
  float ry = 1.0f - fy;
  float rz = 1.0f - fz;
  
  float ryrz = ry*rz;
  float ryfz = ry*fz;
  float fyrz = fy*rz;
  float fyfz = fy*fz;

  do
    {
    vtkResliceRound((rx*(ryrz*inPtr[i000]+ryfz*inPtr[i001]+
			 fyrz*inPtr[i010]+fyfz*inPtr[i011])
		     + fx*(ryrz*inPtr[i100]+ryfz*inPtr[i101]+
			   fyrz*inPtr[i110]+fyfz*inPtr[i111])),
		    *outPtr++);
    inPtr++;
    }
  while (--numscalars);

  return 1;
}			  

// Do nearest-neighbor interpolation of the input data 'inPtr' of extent 
// 'inExt' at the 'point'.  The result is placed at 'outPtr'.  
// If the lookup data is beyond the extent 'inExt', set 'outPtr' to
// the background color 'background'.  
// The number of scalar components in the data is 'numscalars'

template <class T>
static int vtkNearestNeighborInterpolation(float *point, T *inPtr, T *outPtr,
                                           T *background, int numscalars, 
                                           int inExt[6], int inInc[3])
{
  int inIdX = vtkResliceFloor(point[0]+0.5f)-inExt[0];
  int inIdY = vtkResliceFloor(point[1]+0.5f)-inExt[2];
  int inIdZ = vtkResliceFloor(point[2]+0.5f)-inExt[4];

  if (inIdX < 0 || inIdX > inExt[1]-inExt[0]
      || inIdY < 0 || inIdY > inExt[3]-inExt[2]
      || inIdZ < 0 || inIdZ > inExt[5]-inExt[4] )
    {
    if (background)
      {
      vtkCopyPixel(outPtr,background,numscalars);
      }
    return 0;
    }
  else 
    {
    inPtr += inIdX*inInc[0]+inIdY*inInc[1]+inIdZ*inInc[2];
    vtkCopyPixel(outPtr,inPtr,numscalars);

    return 1;
    }
} 

// nearest-neighbor interpolation with wrap-around behaviour
template <class T>
static int vtkNearestNeighborInterpolationRepeat(float *point, T *inPtr, 
						 T *outPtr,
						 T *mirror, int numscalars, 
						 int inExt[6], int inInc[3])
{
  int inIdX = vtkResliceFloor(point[0]+0.5f)-inExt[0];
  int inIdY = vtkResliceFloor(point[1]+0.5f)-inExt[2];
  int inIdZ = vtkResliceFloor(point[2]+0.5f)-inExt[4];

  int inExtX = inExt[1]-inExt[0]+1;
  int inExtY = inExt[3]-inExt[2]+1;
  int inExtZ = inExt[5]-inExt[4]+1;

  if (mirror)
    {
    inIdX = vtkInterpolateMirror(inIdX,inExtX);
    inIdY = vtkInterpolateMirror(inIdY,inExtY);
    inIdZ = vtkInterpolateMirror(inIdZ,inExtZ);
    }
  else
    {
    inIdX = vtkInterpolateWrap(inIdX,inExtX);
    inIdY = vtkInterpolateWrap(inIdY,inExtY);
    inIdZ = vtkInterpolateWrap(inIdZ,inExtZ);
    }
  
  inPtr += inIdX*inInc[0]+inIdY*inInc[1]+inIdZ*inInc[2];
  vtkCopyPixel(outPtr,inPtr,numscalars);

  return 1; 
} 

// Do tricubic interpolation of the input data 'inPtr' of extent 'inExt' 
// at the 'point'.  The result is placed at 'outPtr'.  
// The number of scalar components in the data is 'numscalars'

// The tricubic interpolation ensures that both the intensity and
// the first derivative of the intensity are smooth across the
// image.  The first derivative is estimated using a 
// centered-difference calculation.


// helper function: set up the lookup indices and the interpolation 
// coefficients

void vtkImageResliceSetInterpCoeffs(float F[4],int *l, int *m, float f, 
		     int interpMode)
{   
  float fp1,fm1,fm2;

  switch (interpMode)
    {
    case 7:     // cubic interpolation
      *l = 0; *m = 4; 
      fm1 = f-1;
      F[0] = -f*fm1*fm1/2;
      F[1] = ((3*f-2)*f-2)*fm1/2;
      F[2] = -((3*f-4)*f-1)*f/2;
      F[3] = f*f*fm1/2;
      break;
    case 0:     // no interpolation
    case 2:
    case 4:
    case 6:
      *l = 1; *m = 2; 
      F[1] = 1;
      F[0] = F[2] = F[3] = 0.0f;
      break;
    case 1:     // linear interpolation
      *l = 1; *m = 3;
      F[0] = F[3] = 0.0;
      F[1] = 1-f;
      F[2] = f;
      break;
    case 3:     // quadratic interpolation
      *l = 1; *m = 4; 
      fm1 = f-1; fm2 = fm1-1;
      F[0] = 0.0f;
      F[1] = fm1*fm2/2;
      F[2] = -f*fm2;
      F[3] = f*fm1/2;
      break;
    case 5:     // quadratic interpolation
      *l = 0; *m = 3; 
      fp1 = f+1; fm1 = f-1; 
      F[0] = f*fm1/2;
      F[1] = -fp1*fm1;
      F[2] = fp1*f/2;
      F[3] = 0.0f;
      break;
    }
}

// tricubic interpolation
template <class T>
static int vtkTricubicInterpolation(float *point, T *inPtr, T *outPtr,
				    T *background, int numscalars, 
				    int inExt[6], int inInc[3])
{
  float fx,fy,fz;
  int floorX = vtkResliceFloor(point[0],fx);
  int floorY = vtkResliceFloor(point[1],fy);
  int floorZ = vtkResliceFloor(point[2],fz);

  int inIdX = floorX-inExt[0];
  int inIdY = floorY-inExt[2];
  int inIdZ = floorZ-inExt[4];

  // the doInterpX,Y,Z variables are 0 if interpolation
  // does not have to be done in the specified direction,
  // i.e. if the x, y or z lookup indices have no fractional
  // component.   
  int doInterpX = (fx != 0);
  int doInterpY = (fy != 0);
  int doInterpZ = (fz != 0);

  // check whether we can do cubic interpolation, quadratic, linear, or none
  // in each of the three directions
  if (inIdX < 0 || inIdX+doInterpX > inExt[1]-inExt[0] ||
      inIdY < 0 || inIdY+doInterpY > inExt[3]-inExt[2] ||
      inIdZ < 0 || inIdZ+doInterpZ > inExt[5]-inExt[4])
    {// out of bounds: clear to background color
    if (background)
      {
      vtkCopyPixel(outPtr,background,numscalars);
      }
    return 0;
    }
  else 
    {// do tricubic interpolation
    float fX[4],fY[4],fZ[4];
    float vY,vZ,val;
    T *inPtr1, *inPtr2;
    int i,j,k,l,jl,jm,kl,km,ll,lm;
    int factX[4],factY[4],factZ[4];
    
    // depending on whether we are at the edge of the 
    // input extent, choose the appropriate interpolation
    // method to use

    int interpModeX = ((inIdX > 0) << 2) + 
                      ((inIdX+2 <= inExt[1]-inExt[0]) << 1) +
                      doInterpX;
    int interpModeY = ((inIdY > 0) << 2) + 
                      ((inIdY+2 <= inExt[3]-inExt[2]) << 1) +
                      doInterpY;
    int interpModeZ = ((inIdZ > 0) << 2) + 
	              ((inIdZ+2 <= inExt[5]-inExt[4]) << 1) +
		      doInterpZ;

    vtkImageResliceSetInterpCoeffs(fX,&ll,&lm,fx,interpModeX);
    vtkImageResliceSetInterpCoeffs(fY,&kl,&km,fy,interpModeY);
    vtkImageResliceSetInterpCoeffs(fZ,&jl,&jm,fz,interpModeZ);

    for (i = 0; i < 4; i++)
      {
      factX[i] = (inIdX+i-1)*inInc[0];
      factY[i] = (inIdY+i-1)*inInc[1];
      factZ[i] = (inIdZ+i-1)*inInc[2];
      }

    // set things up so that we can unroll the inner X loop safely
    for (l = 0; l < ll; l++)
      {
      factX[l] = inIdX*inInc[0];
      }
    for (l = lm; l < 4; l++)
      {
      factX[l] = inIdX*inInc[0];
      }

    // Finally, here is the tricubic interpolation
    // (or cubic-cubic-linear, or cubic-nearest-cubic, etc)
    do
      {
      val = 0;
      for (j = jl; j < jm; j++)
	{
	inPtr1 = inPtr + factZ[j];
	vZ = 0;
	for (k = kl; k < km; k++)
	  {
	  inPtr2 = inPtr1 + factY[k];
	  vY = *(inPtr2+factX[0]) * fX[0] +
	       *(inPtr2+factX[1]) * fX[1] +
	       *(inPtr2+factX[2]) * fX[2] +
	       *(inPtr2+factX[3]) * fX[3];
	  vZ += vY*fY[k]; 
	  }
	val += vZ*fZ[j];
	}
      vtkResliceClamp(val,*outPtr++); // clamp to limits of type
      inPtr++;
      }
    while (--numscalars);

    return 1;
    }
}		  

// tricubic interpolation with wrap-around behaviour
template <class T>
static int vtkTricubicInterpolationRepeat(float *point, T *inPtr, T *outPtr,
					  T *mirror, int numscalars, 
					  int inExt[6], int inInc[3])
{
  int factX[4],factY[4],factZ[4];

  float fx,fy,fz;
  int floorX = vtkResliceFloor(point[0],fx);
  int floorY = vtkResliceFloor(point[1],fy);
  int floorZ = vtkResliceFloor(point[2],fz);

  float fX[4],fY[4],fZ[4];
  float vY,vZ,val;
  T *inPtr1, *inPtr2;
  int i,j,k,jl,jm,kl,km;

  int inIdX = floorX-inExt[0];
  int inIdY = floorY-inExt[2];
  int inIdZ = floorZ-inExt[4];

  int inExtX = inExt[1]-inExt[0]+1;
  int inExtY = inExt[3]-inExt[2]+1;
  int inExtZ = inExt[5]-inExt[4]+1;

  if (mirror)
    {
    for (i = 0; i < 4; i++)
      {
      factX[i] = vtkInterpolateMirror(inIdX-1+i,inExtX)*inInc[0];
      factY[i] = vtkInterpolateMirror(inIdY-1+i,inExtY)*inInc[1];
      factZ[i] = vtkInterpolateMirror(inIdZ-1+i,inExtZ)*inInc[2];
      }
    }
  else
    {
    for (i = 0; i < 4; i++)
      {
      factX[i] = vtkInterpolateWrap(inIdX-1+i,inExtX)*inInc[0];
      factY[i] = vtkInterpolateWrap(inIdY-1+i,inExtY)*inInc[1];
      factZ[i] = vtkInterpolateWrap(inIdZ-1+i,inExtZ)*inInc[2];
      }
    }

  vtkImageResliceSetInterpCoeffs(fX,&i,&i,fx,7);
  vtkImageResliceSetInterpCoeffs(fY,&kl,&km,fy,6+(fy != 0));
  vtkImageResliceSetInterpCoeffs(fZ,&jl,&jm,fz,6+(fz != 0));

  // Finally, here is the tricubic interpolation
  do
    {
    val = 0;
    for (j = jl; j < jm; j++)
      {
      inPtr1 = inPtr + factZ[j];
      vZ = 0;
      for (k = kl; k < km; k++)
	{
	inPtr2 = inPtr1 + factY[k];
	vY = *(inPtr2+factX[0]) * fX[0] +
	     *(inPtr2+factX[1]) * fX[1] +
	     *(inPtr2+factX[2]) * fX[2] +
	     *(inPtr2+factX[3]) * fX[3];
	vZ += vY*fY[k]; 
	}
      val += vZ*fZ[j];
      }
    vtkResliceClamp(val,*outPtr++); // clamp to limits of type
    inPtr++;
    }
  while (--numscalars);

  return 1;
}		  

//----------------------------------------------------------------------------
// Some helper functions
//----------------------------------------------------------------------------

// Convert background color from float to appropriate type, or set up
// the pointer to distinguish between Wrap and Mirror

template <class T>
static void vtkAllocBackground(vtkImageReslice *self, T **background_ptr, 
			       int numComponents)
{
  if (self->GetWrap() || self->GetMirror())
    {
    // kludge to differentiate between wrap and mirror
    *background_ptr = (T *)self->GetMirror();
    }
  else
    {
    int i;
    *background_ptr = new T[numComponents];
    T *background = *background_ptr;

    for (i = 0; i < numComponents; i++)
      {
      if (i < 4)
	{
	vtkResliceClamp(self->GetBackgroundColor()[i],background[i]);
	}
      else
	{
	background[i] = 0;
	}
      }
    }
}

template <class T>
static void vtkFreeBackground(vtkImageReslice *self, T **background_ptr)
{
  if (!(self->GetWrap() || self->GetMirror()))
    {
    delete [] *background_ptr;
    }
  *background_ptr = NULL;
}

// get appropriate interpolation function
template <class T>
static void vtkGetResliceInterpFunc(vtkImageReslice *self, 
				    int (**interpolate)(float *point, 
							T *inPtr, T *outPtr,
							T *background, 
							int numscalars, 
							int inExt[6], 
							int inInc[3]))
{
  if (self->GetWrap() || self->GetMirror())
    {
    switch (self->GetInterpolationMode())
      {
      case VTK_RESLICE_NEAREST:
	*interpolate = &vtkNearestNeighborInterpolationRepeat;
	break;
      case VTK_RESLICE_LINEAR:
	*interpolate = &vtkTrilinearInterpolationRepeat;
	break;
      case VTK_RESLICE_CUBIC:
	*interpolate = &vtkTricubicInterpolationRepeat;
	break;
      }
    }
  else
    {
    switch (self->GetInterpolationMode())
      {
      case VTK_RESLICE_NEAREST:
	*interpolate = &vtkNearestNeighborInterpolation;
	break;
      case VTK_RESLICE_LINEAR:
	*interpolate = &vtkTrilinearInterpolation;
	break;
      case VTK_RESLICE_CUBIC:
	*interpolate = &vtkTricubicInterpolation;
	break;
      }
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
  int numscalars;
  int idX, idY, idZ;
  int outIncX, outIncY, outIncZ;
  int inExt[6], inInc[3];
  unsigned long count = 0;
  unsigned long target;
  float point[4];
  float f;
  float *inSpacing,*inOrigin,*outSpacing,*outOrigin,inInvSpacing[3];
  T *background;
  int (*interpolate)(float *point, T *inPtr, T *outPtr,
                     T *background, int numscalars,
                     int inExt[6], int inInc[3]);

  vtkAbstractTransform *transform = self->GetResliceTransform();
  vtkMatrix4x4 *matrix = self->GetResliceAxes();

  inOrigin = inData->GetOrigin();
  inSpacing = inData->GetSpacing();
  outOrigin = self->GetOutputOrigin();
  outSpacing = self->GetOutputSpacing();

  // save effor later: invert inSpacing
  inInvSpacing[0] = 1.0f/inSpacing[0];
  inInvSpacing[1] = 1.0f/inSpacing[1];
  inInvSpacing[2] = 1.0f/inSpacing[2];

  // find maximum input range
  inData->GetExtent(inExt);
  
  target = (unsigned long)
    ((outExt[5]-outExt[4]+1)*(outExt[3]-outExt[2]+1)/50.0);
  target++;
  
  // Get Increments to march through data 
  inData->GetIncrements(inInc);
  outData->GetContinuousIncrements(outExt, outIncX, outIncY, outIncZ);
  numscalars = inData->GetNumberOfScalarComponents();
  
  // set color for area outside of input volume extent
  vtkAllocBackground(self,&background,numscalars);

  // Set interpolation method
  vtkGetResliceInterpFunc(self,&interpolate);

  // Loop through output pixels
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
	point[0] = idX*outSpacing[0] + outOrigin[0];
	point[1] = idY*outSpacing[1] + outOrigin[1];
	point[2] = idZ*outSpacing[2] + outOrigin[2];

	if (matrix)
	  {
	  point[3] = 1.0f;
	  matrix->MultiplyPoint(point,point);
	  f = 1.0f/point[3];
	  point[0] *= f; // deal with w if the matrix
	  point[1] *= f; //   was a Homogeneous transform
	  point[2] *= f;
	  }
	if (transform)
	  {
	  transform->InternalTransformPoint(point,point);
	  }
	  
	point[0] = (point[0] - inOrigin[0])*inInvSpacing[0];
	point[1] = (point[1] - inOrigin[1])*inInvSpacing[1];
	point[2] = (point[2] - inOrigin[2])*inInvSpacing[2];
	
	interpolate(point, inPtr, outPtr, background, 
		    numscalars, inExt, inInc);

	outPtr += numscalars; 
	}
      outPtr += outIncY;
      }
    outPtr += outIncZ;
    }

  vtkFreeBackground(self,&background);
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
  if (this->Optimization && 
      !(this->ResliceTransform && 
	!this->ResliceTransform->IsA("vtkHomogeneousTransform")))
    {
    this->OptimizedThreadedExecute(inData,outData,outExt,id);
    return;
    }

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
    }
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
// The remainder of this file is the 'optimized' version of the code.
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

void vtkImageReslice::OptimizedComputeInputUpdateExtent(int inExt[6], 
							int outExt[6])
{
  int i,j,k;
  int idX,idY,idZ;
  float xAxis[4], yAxis[4], zAxis[4], origin[4];
  float point[4],f;

  int wrap = (this->GetWrap() || this->GetInterpolationMode() != VTK_RESLICE_NEAREST);

  // convert matrix from world coordinates to pixel indices
  vtkMatrix4x4 *matrix = this->GetIndexMatrix();
  for (i = 0; i < 4; i++)
    {
    xAxis[i] = matrix->GetElement(i,0);
    yAxis[i] = matrix->GetElement(i,1);
    zAxis[i] = matrix->GetElement(i,2);
    origin[i] = matrix->GetElement(i,3);
    }

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
      point[j] = point[j] + idX*xAxis[j];
      }
    
    f = 1.0f/point[3];
    point[0] *= f;
    point[1] *= f;
    point[2] *= f;
    
    if (this->GetInterpolationMode() != VTK_RESLICE_NEAREST)
      {
      int extra = (this->GetInterpolationMode() == VTK_RESLICE_CUBIC); 
      for (j = 0; j < 3; j++) 
	{
	k = vtkResliceFloor(point[j])-extra;
	if (k < inExt[2*j])
	  { 
	  inExt[2*j] = k;
	  }
	if (wrap)
	  {
	  k = vtkResliceFloor(point[j])+1+extra;
	  }
	else
	  {
	  k = vtkResliceCeil(point[j])+extra;
	  }
	if (k > inExt[2*j+1])
	  { 
	  inExt[2*j+1] = k;
	  }
	}
      }
    else
      {
      for (j = 0; j < 3; j++) 
	{
	k = vtkResliceFloor(point[j] + 0.5f);
	if (k < inExt[2*j]) 
	  {
	  inExt[2*j] = k;
	  } 
	if (k > inExt[2*j+1])
	  { 
	  inExt[2*j+1] = k;
	  }
	}
      }
    }

  int *wholeExtent = this->GetInput()->GetWholeExtent();
  // Clip, just to make sure we hit _some_ of the input extent
  for (i = 0; i < 3; i++)
    {
    if (inExt[2*i] < wholeExtent[2*i])
      {
      inExt[2*i] = wholeExtent[2*i];
      if (wrap)
	{
	inExt[2*i+1] = wholeExtent[2*i+1];
	}
      }
    if (inExt[2*i+1] > wholeExtent[2*i+1])
      {
      inExt[2*i+1] = wholeExtent[2*i+1];
      if (wrap)
	{
	inExt[2*i] = wholeExtent[2*i];
	}
      }
    if (inExt[2*i] > wholeExtent[2*i+1])
      {
      inExt[2*i] = wholeExtent[2*i+1];
      }
    if (inExt[2*i+1] < wholeExtent[2*i])
      {
      inExt[2*i+1] = wholeExtent[2*i];
      }
    }
}

//----------------------------------------------------------------------------
// helper functions for vtkOptimizedExecute()

// find approximate intersection of line with the plane x = x_min,
// y = y_min, or z = z_min (lower limit of data extent) 

static int intersectionLow(float *point, float *axis, int *sign,
			   int *limit, int ai, int *outExt)
{
  // approximate value of r
  int r;
  float f,p;
  float rd = (limit[ai]*point[3]-point[ai])
    /(axis[ai]-limit[ai]*axis[3]) + 0.5f;
   
  if (rd < outExt[2*ai]) 
    {
    r = outExt[2*ai];
    }
  else if (rd > outExt[2*ai+1])
    {
    r = outExt[2*ai+1];
    }
  else
    {
    r = int(rd);
    }
  
  // move back and forth to find the point just inside the extent
  for (;;)
    {
    f = point[3]+r*axis[3];
    p = point[ai]+r*axis[ai];
    f = 1.0f/f;
    p *= f;
    
    if (vtkResliceFloor(p + 0.5f) < limit[ai])
      {
      r += sign[ai];
      }
    else
      {
      break;
      }
    }

  for (;;)
    {
    f = point[3]+(r-sign[3])*axis[3];
    p = point[ai]+(r-sign[ai])*axis[ai];
    f = 1.0f/f;
    p *= f;
    
    if (vtkResliceFloor(p + 0.5f) >= limit[ai])
      {
      r -= sign[ai];
      }
    else
      {
      break;
      }
    }

  return r;
}

// same as above, but for x = x_max
static int intersectionHigh(float *point, float *axis, int *sign, 
			    int *limit, int ai, int *outExt)
{
  int r;
  float f,p;
  float rd = (limit[ai]*point[3]-point[ai])
      /(axis[ai]-limit[ai]*axis[3]) + 0.5f; 
    
  if (rd < outExt[2*ai])
    { 
    r = outExt[2*ai];
    }
  else if (rd > outExt[2*ai+1])
    {
    r = outExt[2*ai+1];
    }
  else
    {
    r = int(rd);
    }
  
  // move back and forth to find the point just inside the extent
  for (;;)
    {
    f = point[3]+r*axis[3];
    p = point[ai]+r*axis[ai];
    f = 1.0f/f;
    p *= f;
    
    if (vtkResliceFloor(p + 0.5f) > limit[ai])
      {
      r -= sign[ai];
      }
    else
      {
      break;
      }
    }

  for (;;)
    {
    f = point[3]+(r+sign[3])*axis[3];
    p = point[ai]+(r+sign[ai])*axis[ai];
    f = 1.0f/f;
    p *= f;
    
    if (vtkResliceFloor(p + 0.5f) <= limit[ai])
      {
      r += sign[ai];
      }
    else
      {
      break;
      }
    }

  return r;
}

static int isBounded(float *point, float *xAxis, int *inMin, 
		     int *inMax, int ai, int r)
{
  int bi = ai+1; 
  int ci = ai+2;
  if (bi > 2) 
    { 
    bi -= 3; // coordinate index must be 0, 1 or 2 
    } 
  if (ci > 2)
    { 
    ci -= 3;
    }
  float f = point[3]+r*xAxis[3];
  f = 1.0f/f;
  float fbp = point[bi]+r*xAxis[bi];
  float fcp = point[ci]+r*xAxis[ci];
  fbp *= f;
  fcp *= f;

  int bp = vtkResliceFloor(fbp + 0.5f);
  int cp = vtkResliceFloor(fcp + 0.5f);
  
  return (bp >= inMin[bi] && bp <= inMax[bi] &&
	  cp >= inMin[ci] && cp <= inMax[ci]);
}

// this huge mess finds out where the current output raster
// line intersects the input volume 
int vtkImageReslice::FindExtent(int& r1, int& r2, float *point, float *xAxis, 
				int *inMin, int *inMax, int *outExt)
{
  int i, ix, iy, iz;
  int sign[3];
  int indx1[4],indx2[4];
  float f1,f2,p1,p2;

  // find signs of components of x axis 
  // (this is complicated due to the homogeneous coordinate)
  for (i = 0; i < 3; i++)
    {
    f1 = point[3];
    p1 = point[i];
    f1 = 1.0f/f1;
    p1 *= f1;

    f2 = point[3]+xAxis[3];
    p2 = point[i]+xAxis[i];
    f2 = 1.0f/f2;
    p2 *= f2;

    if (p1 <= p2)
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

  r1 = intersectionLow(point,xAxis,sign,inMin,ix,outExt);
  r2 = intersectionHigh(point,xAxis,sign,inMax,ix,outExt);
  
  // find points of intersections
  // first, find w-value for perspective (will usually be 1)
  f1 = point[3]+r1*xAxis[3];
  f2 = point[3]+r2*xAxis[3];
  f1 = 1.0f/f1;
  f2 = 1.0f/f2;
  
  for (i = 0; i < 3; i++)
    {
    p1 = point[i]+r1*xAxis[i];
    p2 = point[i]+r2*xAxis[i];
    p1 *= f1;
    p2 *= f2;

    indx1[i] = vtkResliceFloor(p1 + 0.5f);
    indx2[i] = vtkResliceFloor(p2 + 0.5f);
    }
  if (isBounded(point,xAxis,inMin,inMax,ix,r1))
    { // passed through x face, check opposing face
    if (isBounded(point,xAxis,inMin,inMax,ix,r2))
      {
      return sign[ix];
      }
    
    if (indx2[iy] < inMin[iy])
      { // check y face
      r2 = intersectionLow(point,xAxis,sign,inMin,iy,outExt);
      if (isBounded(point,xAxis,inMin,inMax,iy,r2))
	{
	return sign[ix];
	}
      }
    else if (indx2[iy] > inMax[iy])
      { // check other y face
      r2 = intersectionHigh(point,xAxis,sign,inMax,iy,outExt);
      if (isBounded(point,xAxis,inMin,inMax,iy,r2))
	{
	return sign[ix];
	}
      }
    
    if (indx2[iz] < inMin[iz])
      { // check z face
      r2 = intersectionLow(point,xAxis,sign,inMin,iz,outExt);
      if (isBounded(point,xAxis,inMin,inMax,iz,r2))
	{
	return sign[ix];
	}
      }
    else if (indx2[iz] > inMax[iz])
      { // check other z face
      r2 = intersectionHigh(point,xAxis,sign,inMax,iz,outExt);
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
	r1 = intersectionLow(point,xAxis,sign,inMin,iy,outExt);
	if (isBounded(point,xAxis,inMin,inMax,iy,r1))
	  {
	  return sign[ix];
	  }
	}
    else if (indx1[iy] > inMax[iy])
      { // check other y face
      r1 = intersectionHigh(point,xAxis,sign,inMax,iy,outExt);
      if (isBounded(point,xAxis,inMin,inMax,iy,r1))
	{
	return sign[ix];
	}
      }
    
    if (indx1[iz] < inMin[iz])
      { // check z face
      r1 = intersectionLow(point,xAxis,sign,inMin,iz,outExt);
      if (isBounded(point,xAxis,inMin,inMax,iz,r1))
	{
	return sign[ix];
	}
      }
    else if (indx1[iz] > inMax[iz])
      { // check other z face
      r1 = intersectionHigh(point,xAxis,sign,inMax,iz,outExt);
      if (isBounded(point,xAxis,inMin,inMax,iz,r1))
	{
	return sign[ix];
	}
      }
    }
  
  if ((indx1[iy] >= inMin[iy] && indx2[iy] < inMin[iy]) ||
      (indx1[iy] < inMin[iy] && indx2[iy] >= inMin[iy]))
    { // line might pass through bottom face
    r1 = intersectionLow(point,xAxis,sign,inMin,iy,outExt);
    if (isBounded(point,xAxis,inMin,inMax,iy,r1))
      {
      if ((indx1[iy] <= inMax[iy] && indx2[iy] > inMax[iy]) ||
	  (indx1[iy] > inMax[iy] && indx2[iy] <= inMax[iy]))
	{ // line might pass through top face
	r2 = intersectionHigh(point,xAxis,sign,inMax,iy,outExt);
	if (isBounded(point,xAxis,inMin,inMax,iy,r2))
	  {
	  return sign[iy];
	  }
	}
      
      if (indx1[iz] < inMin[iz] && indx2[iy] < inMin[iy] ||
	  indx2[iz] < inMin[iz] && indx1[iy] < inMin[iy])
	{ // line might pass through in-to-screen face
	r2 = intersectionLow(point,xAxis,sign,inMin,iz,outExt);
	if (isBounded(point,xAxis,inMin,inMax,iz,r2))
	  {
	  return sign[iy];
	  }
	}
      else if (indx1[iz] > inMax[iz] && indx2[iy] < inMin[iy] ||
	       indx2[iz] > inMax[iz] && indx1[iy] < inMin[iy])
	{ // line might pass through out-of-screen face
	r2 = intersectionHigh(point,xAxis,sign,inMax,iz,outExt);
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
    r2 = intersectionHigh(point,xAxis,sign,inMax,iy,outExt);
    if (isBounded(point,xAxis,inMin,inMax,iy,r2))
      {
      if (indx1[iz] < inMin[iz] && indx2[iy] > inMax[iy] ||
	  indx2[iz] < inMin[iz] && indx1[iy] > inMax[iy])
	{ // line might pass through in-to-screen face
	r1 = intersectionLow(point,xAxis,sign,inMin,iz,outExt);
	if (isBounded(point,xAxis,inMin,inMax,iz,r1))
	  {
	  return sign[iy];
	  }
	}
      else if (indx1[iz] > inMax[iz] && indx2[iy] > inMax[iy] || 
	       indx2[iz] > inMax[iz] && indx1[iy] > inMax[iy])
	{ // line might pass through out-of-screen face
	r1 = intersectionHigh(point,xAxis,sign,inMax,iz,outExt);
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
    r1 = intersectionLow(point,xAxis,sign,inMin,iz,outExt);
    if (isBounded(point,xAxis,inMin,inMax,iz,r1))
      {
      if (indx1[iz] > inMax[iz] || indx2[iz] > inMax[iz])
	{ // line might pass through out-of-screen face
	r2 = intersectionHigh(point,xAxis,sign,inMax,iz,outExt);
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
				int outExt[6], int id, vtkMatrix4x4 *matrix)
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
  int r1,r2;
  float inPoint0[4];
  float inPoint1[4];
  float xAxis[4], yAxis[4], zAxis[4], origin[4];
  float inPoint[4],f;
  T *background;
  int (*interpolate)(float *point, T *inPtr, T *outPtr,
                     T *background, int numscalars,
                     int inExt[6], int inInc[3]);

  // find maximum input range
  self->GetInput()->GetExtent(inExt);

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
  
  // break matrix into a set of axes plus an origin
  // (this allows us to calculate the transform Incrementally)
  for (i = 0; i < 4; i++)
    {
    xAxis[i] = matrix->GetElement(i,0);
    yAxis[i] = matrix->GetElement(i,1);
    zAxis[i] = matrix->GetElement(i,2);
    origin[i] = matrix->GetElement(i,3);
    }

  // set color for area outside of input volume extent
  vtkAllocBackground(self,&background,numscalars);

  // Set interpolation method
  vtkGetResliceInterpFunc(self,&interpolate);

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
      
      if (self->GetWrap() || self->GetMirror())
	{ // wrap-pad like behaviour
	for (idX = outExt[0]; idX <= outExt[1]; idX++)
	  {
	  inPoint[0] = inPoint1[0]+idX*xAxis[0];
	  inPoint[1] = inPoint1[1]+idX*xAxis[1];
	  inPoint[2] = inPoint1[2]+idX*xAxis[2];
	  if (inPoint1[3] != 1.0f || xAxis[3] != 0.0f)
	    { // only do perspective if necessary
	    inPoint[3] = inPoint1[3]+idX*xAxis[3];

	    f = 1.0f/inPoint[3];
	    inPoint[0] *= f;
	    inPoint[1] *= f;
	    inPoint[2] *= f;
	    }

	  interpolate(inPoint, inPtr, outPtr, background, 
		      numscalars, inExt, inInc);
	  outPtr += numscalars;
	  }
	}
      else
	{
	// find intersections of x raster line with the input extent
	if (self->FindExtent(r1,r2,inPoint1,xAxis,inMin,inMax,outExt) < 0)
	  {
	  i = r1;
	  r1 = r2;
	  r2 = i;
	  }

	// bound r1,r2 within reasonable limits
	if (r1 < outExt[0]) 
	  {
	  r1 = outExt[0];
	  }
	if (r2 > outExt[1]) 
	  {
	  r2 = outExt[1];
	  }
	if (r1 > r2) 
	  {
	  r1 = outExt[0];
	  r2 = outExt[0]-1;
	  }

	// clear pixels to left of input extent
	if (numscalars == 1) // optimize for single scalar
	  {
	  for (idX = outExt[0]; idX < r1; idX++) 
	    {
	    *outPtr++ = background[0];
	    }
	  }
	else             // multiple scalars
	  {
	  for (idX = outExt[0]; idX < r1; idX++)
	    {
	    vtkCopyPixel(outPtr,background,numscalars);
	    }
	  }
	
	if (self->GetInterpolationMode() != VTK_RESLICE_NEAREST)
	  { // Trilinear or tricubic
	  for (idX = r1; idX <= r2; idX++)
	    {
	    inPoint[0] = inPoint1[0]+idX*xAxis[0];
	    inPoint[1] = inPoint1[1]+idX*xAxis[1];
	    inPoint[2] = inPoint1[2]+idX*xAxis[2];
	    if (inPoint1[3] != 1.0f || xAxis[3] != 0.0f)
	      { // only do perspective if necessary
	      inPoint[3] = inPoint1[3]+idX*xAxis[3];

	      f = 1.0f/inPoint[3];
	      inPoint[0] *= f;
	      inPoint[1] *= f;
	      inPoint[2] *= f;
	      }
	    
	    interpolate(inPoint, inPtr, outPtr, background, 
			numscalars, inExt, inInc);
	    outPtr += numscalars;
	    }
	  }
	else if (inPoint1[3] != 1.0f || xAxis[3] != 0.0f)
	  {  // Nearest-Neighbor, no extent checks, perspective
	  T *inPtr1;

	  for (idX = r1; idX <= r2; idX++)
	    {
	    inPoint[0] = inPoint1[0]+idX*xAxis[0];
	    inPoint[1] = inPoint1[1]+idX*xAxis[1];
	    inPoint[2] = inPoint1[2]+idX*xAxis[2];
	    inPoint[3] = inPoint1[3]+idX*xAxis[3];
	    
	    f = 1.0f/inPoint[3];
	    inPoint[0] *= f;
	    inPoint[1] *= f;
	    inPoint[2] *= f;
	    
	    inIdX = vtkResliceFloor(inPoint[0]+0.5f)-inExt[0];
	    inIdY = vtkResliceFloor(inPoint[1]+0.5f)-inExt[2];
	    inIdZ = vtkResliceFloor(inPoint[2]+0.5f)-inExt[4];

	    inPtr1 = inPtr+inIdX*inInc[0]+inIdY*inInc[1]+inIdZ*inInc[2];
	    vtkCopyPixel(outPtr,inPtr1,numscalars);
	    }
	  }
	else
	  { // Nearest-neighbor, no extent checks, linear
	  T *inPtr1;

	  for (idX = r1; idX <= r2; idX++)
	    {
	    inPoint[0] = inPoint1[0]+idX*xAxis[0];
	    inPoint[1] = inPoint1[1]+idX*xAxis[1];
	    inPoint[2] = inPoint1[2]+idX*xAxis[2];

	    inIdX = vtkResliceFloor(inPoint[0]+0.5f)-inExt[0];
	    inIdY = vtkResliceFloor(inPoint[1]+0.5f)-inExt[2];
	    inIdZ = vtkResliceFloor(inPoint[2]+0.5f)-inExt[4];

	    inPtr1 = inPtr+inIdX*inInc[0]+inIdY*inInc[1]+inIdZ*inInc[2];
	    vtkCopyPixel(outPtr,inPtr1,numscalars);
	    }
	  }
  
	// clear pixels to right of input extent
	if (numscalars == 1) // optimize for single scalar
	  {
	  for (idX = r2+1; idX <= outExt[1]; idX++)
	    { 
	    *outPtr++ = background[0];
	    }
	  }
	else // multiple scalars
	  {
	  for (idX = r2+1; idX <= outExt[1]; idX++)
	    {
	    vtkCopyPixel(outPtr,background,numscalars);
	    }
	  }
	}

      outPtr += outIncY;
      }
    outPtr += outIncZ;
    }
  
  vtkFreeBackground(self,&background);
}

// vtkOptimizedPermuteExecute is specifically optimized for
// cases where the IndexMatrix has only one non-zero component
// per row, i.e. when the matrix is permutation+scale+translation.
// All of the interpolation coefficients are calculated ahead
// of time instead of on a pixel-by-pixel basis.

template <class T>
static void vtkOptimizedPermuteExecuteLinear(vtkImageReslice *self,
					     vtkImageData *inData, T *inPtr,
					     vtkImageData *outData, T *outPtr,
					     int outExt[6], int id,
					     vtkMatrix4x4 *matrix)
{
  int i, j, k, numscalars;
  int idX, idY, idZ;
  int outIncX, outIncY, outIncZ;
  int inExt[6];
  int inInc[3];
  int clipExt[6];
  unsigned long count = 0;
  unsigned long target;
  int r1,r2;
  T *background;

  // find maximum input range
  self->GetInput()->GetExtent(inExt);

  target = (unsigned long)
    ((outExt[5]-outExt[4]+1)*(outExt[3]-outExt[2]+1)/50.0);
  target++;
  
  // Get Increments to march through data 
  inData->GetIncrements(inInc);
  outData->GetContinuousIncrements(outExt, outIncX, outIncY, outIncZ);
  numscalars = inData->GetNumberOfScalarComponents();

  // set color for area outside of input volume extent
  vtkAllocBackground(self,&background,numscalars);

  for (i = 0; i < 3; i++)
    {
    clipExt[2*i] = 0;
    clipExt[2*i+1] = outExt[2*i+1]-outExt[2*i];
    }

  int *traversal[3];
  float *constants[3];
  float newmat[4][4];
  int region;

  for (j = 0; j < 4; j++)
    {
    for (i = 0; i < 4; i++)
      {
      newmat[i][j] = matrix->GetElement(i,j);
      }
    } 

  int trunc, inId0, inId1, inExtK, outExtJ, doInterp;
  float point,f;
  
  // set up input traversal table for linear interpolation  
  for (j = 0; j < 3; j++)
    {
    outExtJ = outExt[2*j+1]-outExt[2*j]+1;
    traversal[j] = new int[outExtJ*2];
    constants[j] = new float[outExtJ*2];

    for (k = 0; k < 3; k++)
      { // set k to the element which is nonzero
      if (newmat[k][j] != 0)
	{
	break;
	}
      }    
    inExtK = inExt[2*k+1]-inExt[2*k]+1;

    region = 0;
    for (i = 0; i < outExtJ; i++)
      {
      point = newmat[k][3]+(i+outExt[2*j])*newmat[k][j];
      trunc = vtkResliceFloor(point,f);
      constants[j][2*i] = 1.0f-f;
      constants[j][2*i+1] = f;
      
      doInterp = (f != 0);
      inId0 = trunc - inExt[2*k];
      inId1 = inId0+doInterp;

        if (inId0 < 0 || inId1 >= inExtK)
	  {
	  if (region == 1)
	    { // leaving the input extent
	    region = 2;
	    clipExt[2*j+1] = i-1;
	    }
	  }
        else 
	  {
	  if (region == 0)
	    { // entering the input extent
	    region = 1;
	    clipExt[2*j] = i;	   
	    }
	  }
      traversal[j][2*i] = inId0*inInc[k];
      traversal[j][2*i+1] = inId1*inInc[k];
      }
    if (region == 0)
      { // never entered input extent!
      clipExt[2*j] = clipExt[2*j+1]+1;
      }
    }

  int outExtX = outExt[1]-outExt[0]+1;
  int outExtY = outExt[3]-outExt[2]+1;
  int outExtZ = outExt[5]-outExt[4]+1;

  // Loop through output pixels
  for (idZ = 0; idZ < outExtZ; idZ++)
    {
    int idZ0 = 2*idZ;
    int idZ1 = idZ0+1;

    int i0 = traversal[2][idZ0]; 
    int i1 = traversal[2][idZ1];

    float rz = constants[2][idZ0];
    float fz = constants[2][idZ1];
    
    for (idY = 0; idY < outExtY; idY++)
      {
      int idY0 = 2*idY;
      int idY1 = idY0+1;

      int i00 = traversal[1][idY0] + i0;
      int i01 = traversal[1][idY0] + i1;
      int i10 = traversal[1][idY1] + i0;
      int i11 = traversal[1][idY1] + i1;

      float ry = constants[1][idY0];
      float fy = constants[1][idY1];

      float ryrz = ry*rz;
      float ryfz = ry*fz;
      float fyrz = fy*rz;
      float fyfz = fy*fz;

      if (!id) 
	{
	if (!(count%target)) 
	  {
	  self->UpdateProgress(count/(50.0*target));
	  }
	count++;
	}

      // do extent check
      if (idZ < clipExt[4] || idZ > clipExt[5] ||
	  idY < clipExt[2] || idY > clipExt[3])
	{
	r1 = outExtX;
	r2 = outExtX-1;
	}
      else
	{
	r1 = clipExt[0];
	r2 = clipExt[1];
	}
  
      // clear pixels to left of input extent
      for (idX = 0; idX < r1; idX++)
	{
	vtkCopyPixel(outPtr,background,numscalars);
	}

      for (idX = r1; idX <= r2; idX++)
	{
        int idX0 = 2*idX;
	int idX1 = idX0+1;
	
	int t0 = traversal[0][idX0];
	int t1 = traversal[0][idX1];
	
	int i000 = t0 + i00; 
	int i001 = t0 + i01; 
	int i010 = t0 + i10; 
	int i011 = t0 + i11; 
	int i100 = t1 + i00; 
	int i101 = t1 + i01; 
	int i110 = t1 + i10; 
	int i111 = t1 + i11; 
	
	float rx = constants[0][idX0];
	float fx = constants[0][idX1];
	
	T *inPtr0 = inPtr;

	i = numscalars;
	do 
	  {
	  vtkResliceRound((rx*(ryrz*inPtr0[i000]+ryfz*inPtr0[i001]+
			       fyrz*inPtr0[i010]+fyfz*inPtr0[i011])
			   + fx*(ryrz*inPtr0[i100]+ryfz*inPtr0[i101]+
				 fyrz*inPtr0[i110]+fyfz*inPtr0[i111])),
			  *outPtr++);
	  inPtr0++;
	  }
	while (--i);
	}

      // clear pixels to right of input extent
      for (idX = r2+1; idX < outExtX; idX++)
        {
	vtkCopyPixel(outPtr,background,numscalars);
	}

      outPtr += outIncY;
      }
    outPtr += outIncZ;
    }

  for (j = 0; j < 3; j++)
    {
    delete [] traversal[j];
    delete [] constants[j];
    }

  vtkFreeBackground(self,&background);
}

template <class T>
static void vtkOptimizedPermuteExecuteCubic(vtkImageReslice *self,
                                            vtkImageData *inData, T *inPtr,
                                            vtkImageData *outData, T *outPtr,
					    int outExt[6], int id,
					    vtkMatrix4x4 *matrix)
{
  int i, j, k, l, numscalars;
  int idX, idY, idZ;
  int outIncX, outIncY, outIncZ;
  int inExt[6];
  int inInc[3];
  int clipExt[6];
  unsigned long count = 0;
  unsigned long target;
  int r1,r2;
  T *background;

  // find maximum input range
  self->GetInput()->GetExtent(inExt);

  target = (unsigned long)
    ((outExt[5]-outExt[4]+1)*(outExt[3]-outExt[2]+1)/50.0);
  target++;
  
  // Get Increments to march through data 
  inData->GetIncrements(inInc);
  outData->GetContinuousIncrements(outExt, outIncX, outIncY, outIncZ);
  numscalars = inData->GetNumberOfScalarComponents();

  // set color for area outside of input volume extent
  vtkAllocBackground(self,&background,numscalars);

  for (i = 0; i < 3; i++)
    {
    clipExt[2*i] = 0;
    clipExt[2*i+1] = outExt[2*i+1]-outExt[2*i];
    }

  int *traversal[3];
  int *low[3];
  int *high[3];
  float *constants[3];
  float newmat[4][4];
  int region;

  for (j = 0; j < 4; j++)
    {
    for (i = 0; i < 4; i++)
      {
      newmat[i][j] = matrix->GetElement(i,j);
      }
    } 

  int trunc, inId[4], inExtK, outExtJ, doInterp, interpMode;
  float point,f;
  
  // set up input traversal table for cubic interpolation  
  for (j = 0; j < 3; j++)
    {
    outExtJ = outExt[2*j+1]-outExt[2*j]+1;
    traversal[j] = new int[outExtJ*4];
    constants[j] = new float[outExtJ*4];
    low[j] = new int[outExtJ];
    high[j] = new int[outExtJ];

    for (k = 0; k < 3; k++)
      { // set k to the element which is nonzero
      if (newmat[k][j] != 0)
	{
	break;
	}
      }
    inExtK = inExt[2*k+1]-inExt[2*k]+1;

    region = 0;
    for (i = 0; i < outExtJ; i++)
      {
      point = newmat[k][3]+(i+outExt[2*j])*newmat[k][j];
      trunc = vtkResliceFloor(point,f);
      doInterp = (f != 0);
      inId[1] = trunc - inExt[2*k];
      inId[0] = inId[1]-1;
      inId[2] = inId[1]+1;
      inId[3] = inId[1]+2;

      if (self->GetMirror())
	{
	inId[0] = vtkInterpolateMirror(inId[0], inExtK);
	inId[1] = vtkInterpolateMirror(inId[1], inExtK);
	inId[2] = vtkInterpolateMirror(inId[2], inExtK);
	inId[3] = vtkInterpolateMirror(inId[3], inExtK);
	interpMode = 6+doInterp;
	region = 1;
	}
      else if (self->GetWrap())
	{
	inId[0] = vtkInterpolateWrap(inId[0], inExtK);
	inId[1] = vtkInterpolateWrap(inId[1], inExtK);
	inId[2] = vtkInterpolateWrap(inId[2], inExtK);
	inId[3] = vtkInterpolateWrap(inId[3], inExtK);
	interpMode = 6+doInterp;
	region = 1;
	}      
      else
	{
        if (inId[1] < 0 || inId[1]+doInterp >= inExtK)
	  {
	  if (region == 1)
	    { // leaving the input extent
	    region = 2;
	    clipExt[2*j+1] = i-1;
	    }
	  }
        else 
	  {
	  if (region == 0)
	    { // entering the input extent
	    region = 1;
	    clipExt[2*j] = i;	   
	    }
	  }
	interpMode = ((inId[1] > 0) << 2) + 
	              ((inId[3] < inExtK) << 1) +
	              doInterp; 
	}
      vtkImageResliceSetInterpCoeffs(&constants[j][4*i],&low[j][i],
				     &high[j][i],f,interpMode);
      // set default values
      for (l = 0; l < 4; l++)
	{
	traversal[j][4*i+l] = inId[1]*inInc[k];
	}
      for (l = low[j][i]; l < high[j][i]; l++)
	{ 
	traversal[j][4*i+l] = inId[l]*inInc[k];
	}
      }
    if (region == 0)
      { // never entered input extent!
      clipExt[2*j] = clipExt[2*j+1]+1;
      }
    }

  int outExtX = outExt[1]-outExt[0]+1;
  int outExtY = outExt[3]-outExt[2]+1;
  int outExtZ = outExt[5]-outExt[4]+1;

  // Loop through output pixels
  for (idZ = 0; idZ < outExtZ; idZ++)
    {
    int idZ0 = idZ*4;
    int lz = low[2][idZ];
    int hz = high[2][idZ];
    float fZ[4];
    int iZ[4];

    for (i = 0; i < 4; i++)
      {
      fZ[i] = constants[2][idZ0+i];
      iZ[i] = traversal[2][idZ0+i];
      }

    for (idY = 0; idY < outExtY; idY++)
      {
      int idY0 = idY*4;
      int ly = low[1][idY];
      int hy = high[1][idY];
      float fY[4];
      int iY[4];
      
      for (i = 0; i < 4; i++)
	{
        fY[i] = constants[1][idY0+i];
	iY[i] = traversal[1][idY0+i];
	}

      if (!id) 
	{
	if (!(count%target)) 
	  {
	  self->UpdateProgress(count/(50.0*target));
	  }
	count++;
	}

      // do extent check
      if (idZ < clipExt[4] || idZ > clipExt[5] ||
	  idY < clipExt[2] || idY > clipExt[3])
	{
	r1 = outExtX;
	r2 = outExtX-1;
	}
      else
	{
	r1 = clipExt[0];
	r2 = clipExt[1];
	}
  
      // clear pixels to left of input extent
      for (idX = 0; idX < r1; idX++)
	{
	vtkCopyPixel(outPtr,background,numscalars);
	}

      for (idX = r1; idX <= r2; idX++)
	{
	int idX0 = idX*4;
	float fX[4];
	int iX[4];

	for (i = 0; i < 4; i++)
	  {
	  fX[i] = constants[0][idX0+i];
	  iX[i] = traversal[0][idX0+i];
	  }
	
	T *inPtr0 = inPtr;
	T *inPtr1,*inPtr2;
	float val,vY,vZ;

	l = numscalars;
	do
	  {
	  val = 0;
	  for (k = lz; k < hz; k++)
	    {
	    inPtr1 = inPtr0 + iZ[k];
	    vZ = 0;
	    for (j = ly; j < hy; j++)
	      {
	      inPtr2 = inPtr1 + iY[j];
	      vY = *(inPtr2+iX[0]) * fX[0] +
		   *(inPtr2+iX[1]) * fX[1] +
		   *(inPtr2+iX[2]) * fX[2] +
		   *(inPtr2+iX[3]) * fX[3];
	      vZ += vY*fY[j]; 
	      }
	    val += vZ*fZ[k];
	    }
	  vtkResliceClamp(val,*outPtr++); // clamp to limits of type
	  inPtr0++;
	  }
	while (--l);
	}
      
      // clear pixels to right of input extent
      for (idX = r2+1; idX < outExtX; idX++)
        {
	vtkCopyPixel(outPtr,background,numscalars);
	}

      outPtr += outIncY;
      }
    outPtr += outIncZ;
    }

  for (j = 0; j < 3; j++)
    {
    delete [] traversal[j];
    delete [] constants[j];
    delete [] high[j];
    delete [] low[j];
    }

  vtkFreeBackground(self,&background);
}

template <class T>
static void vtkOptimizedPermuteExecuteNearest(vtkImageReslice *self,
					      vtkImageData *inData, T *inPtr,
					      vtkImageData *outData, T *outPtr,
					      int outExt[6], int id,
					      vtkMatrix4x4 *matrix)
{
  int i, j, k, numscalars;
  int idX, idY, idZ;
  int outIncX, outIncY, outIncZ;
  int inExt[6];
  int inInc[3];
  int clipExt[6];
  unsigned long count = 0;
  unsigned long target;
  int r1,r2;
  T *background,*inPtr0,*inPtr1,*inPtr2;

  // find maximum input range
  self->GetInput()->GetExtent(inExt);

  target = (unsigned long)
    ((outExt[5]-outExt[4]+1)*(outExt[3]-outExt[2]+1)/50.0);
  target++;
  
  // Get Increments to march through data 
  inData->GetIncrements(inInc);
  outData->GetContinuousIncrements(outExt, outIncX, outIncY, outIncZ);
  numscalars = inData->GetNumberOfScalarComponents();

  // set color for area outside of input volume extent
  vtkAllocBackground(self,&background,numscalars);

  for (i = 0; i < 3; i++)
    {
    clipExt[2*i] = 0;
    clipExt[2*i+1] = outExt[2*i+1]-outExt[2*i];
    }

  int *traversal[3];
  float newmat[4][4];
  int region;

  for (j = 0; j < 4; j++)
    {
    for (i = 0; i < 4; i++)
      {
      newmat[i][j] = matrix->GetElement(i,j);
      }
    } 

  int inId, inExtK, outExtJ;

  // set up input traversal table for nearest-neighbor interpolation  
  for (j = 0; j < 3; j++)
    {
    outExtJ = outExt[2*j+1]-outExt[2*j]+1;
    traversal[j] = new int[outExtJ];

    for (k = 0; k < 3; k++)
      { // set k to the element which is nonzero
      if (newmat[k][j] != 0)
	{
	break;
	}
      }    
    inExtK = inExt[2*k+1]-inExt[2*k]+1;

    region = 0;
    for (i = 0; i < outExtJ; i++)
      {
      inId = vtkResliceFloor((newmat[k][3]+(i+outExt[2*j])*newmat[k][j])+0.5f)
	           - inExt[2*k];
      if (self->GetMirror())
	{
	inId = vtkInterpolateMirror(inId, inExtK);
	region = 1;
	}
      else if (self->GetWrap())
	{
	inId = vtkInterpolateWrap(inId, inExtK);
	region = 1;
	}
      else
	{
        if (inId < 0 || inId >= inExtK)
	  {
  	  if (region == 1)
	    { // leaving the input extent
            region = 2;
	    clipExt[2*j+1] = i-1;
	    }
 	  }
        else 
	  {
	  if (region == 0)
	    { // entering the input extent
	    region = 1;
	    clipExt[2*j] = i;	   
	    }
	  }
        }
      traversal[j][i] = inId*inInc[k];
      }
    if (region == 0)
      { // never entered input extent!
      clipExt[2*j] = clipExt[2*j+1]+1;
      }
    }
  
  int outExtX = outExt[1]-outExt[0]+1;
  int outExtY = outExt[3]-outExt[2]+1;
  int outExtZ = outExt[5]-outExt[4]+1;

  // Loop through output pixels
  for (idZ = 0; idZ < outExtZ; idZ++)
    {
    inPtr0 = inPtr + traversal[2][idZ]; 
    
    for (idY = 0; idY < outExtY; idY++)
      {
      inPtr1 = inPtr0 + traversal[1][idY];

      if (!id) 
	{
	if (!(count%target)) 
	  {
	  self->UpdateProgress(count/(50.0*target));
	  }
	count++;
	}

      // do extent check
      if (idZ < clipExt[4] || idZ > clipExt[5] ||
	  idY < clipExt[2] || idY > clipExt[3])
	{
	r1 = outExtX;
	r2 = outExtX-1;
	}
      else
	{
	r1 = clipExt[0];
	r2 = clipExt[1];
	}
 
      if (numscalars == 1) // fast path for single-scalar data
	{ 
        // clear pixels to left of input extent
	for (idX = 0; idX < r1; idX++)
	  {
	  *outPtr++ = background[0];
	  }

	// do nearest-neighbor interpolation
	for (idX = r1; idX <= r2; idX++)
	  {
	  *outPtr++ = inPtr1[traversal[0][idX]];
	  }

	// clear pixels to right of input extent
	for (idX = r2+1; idX < outExtX; idX++)
	  {
	  *outPtr++ = background[0];
	  }
	}
      else // multiple scalars
	{
        // clear pixels to left of input extent
        for (idX = outExt[0]; idX < r1; idX++)
 	  {
	  vtkCopyPixel(outPtr,background,numscalars);
	  }

	// do nearest-neighbor interpolation
        for (idX = r1; idX <= r2; idX++)
	  {
	  inPtr2 = inPtr1 + traversal[0][idX];
	  vtkCopyPixel(outPtr,inPtr2,numscalars);
	  }

        // clear pixels to right of input extent
        for (idX = r2+1; idX < outExtX; idX++)
          {
	  vtkCopyPixel(outPtr,background,numscalars);
	  }
	}

      outPtr += outIncY;
      }
    outPtr += outIncZ;
    }

  for (j = 0; j < 3; j++)
    {
    delete [] traversal[j];
    }

  vtkFreeBackground(self,&background);
}

//----------------------------------------------------------------------------
// check a matrix to see whether it is the identity matrix

static int vtkIsIdentityMatrix(vtkMatrix4x4 *matrix)
{
  static double identity[16] = {1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1};
  int i,j;

  for (i = 0; i < 4; i++)
    {
    for (j = 0; j < 4; j++)
      {
      if (matrix->GetElement(i,j) != identity[4*i+j])
	{
	return 0;
	}
      }
    }
  return 1;
}

// check a matrix to ensure that it is a permutation+scale+translation
// matrix

static int vtkIsPermutationMatrix(vtkMatrix4x4 *matrix)
{
  int i,j,k;

  for (i = 0; i < 3; i++)
    {
    if (matrix->GetElement(3,i) != 0.0)
      {
      return 0;
      }
    }
  if (matrix->GetElement(3,3) != 1.0)
    {
    return 0;
    }
  for (j = 0; j < 3; j++)
    {
    k = 0;
    for (i = 0; i < 3; i++)
      {
      if (matrix->GetElement(i,j) != 0.0)
	{
	k++;
	}
      }
    if (k != 1)
      {
      return 0;
      }
    }
  return 1;
}

// Check to see if we can do nearest-neighbor instead of linear or cubic.  
// This check only works on permutation+scale+translation matrices.

static inline int vtkCanUseNearestNeighbor(vtkMatrix4x4 *matrix, int outExt[6])
{
  int i,j;
  double x,y;

  // loop through dimensions
  for (i = 0; i < 3; i++)
    {
    x = 0;
    for (j = 0; j < 3; j++)
      {
      x += matrix->GetElement(i,j);
      }
    y = matrix->GetElement(i,3);
    if (outExt[2*i] == outExt[2*i+1])
      {
      y += x*outExt[2*i];
      x = 0;
      }
    if (x != int(x) || y != int(y)) 
      {
      return 0;
      }
    }
  return 1;
}

// the OptimizedPermuteExecute path is taken when the output slices are
// orthogonal to the input slices

template <class T>
static void vtkOptimizedPermuteExecute(vtkImageReslice *self,
				       vtkImageData *inData, T *inPtr,
				       vtkImageData *outData, T *outPtr,
				       int outExt[6], int id,
				       vtkMatrix4x4 *matrix)
{
  if (self->GetInterpolationMode() == VTK_RESLICE_NEAREST ||
      vtkCanUseNearestNeighbor(matrix,outExt))
    {
    vtkOptimizedPermuteExecuteNearest(self,inData,inPtr,outData,outPtr,
				      outExt,id,matrix);
    }
  else if (self->GetInterpolationMode() == VTK_RESLICE_LINEAR)
    {
    vtkOptimizedPermuteExecuteLinear(self,inData,inPtr,outData,outPtr,
				     outExt,id,matrix);
    }
  else if (self->GetInterpolationMode() == VTK_RESLICE_CUBIC)
    {
    vtkOptimizedPermuteExecuteCubic(self,inData,inPtr,outData,outPtr,
				    outExt,id,matrix);
    }
}
    
//----------------------------------------------------------------------------
// The transform matrix supplied by the user converts output coordinates
// to input coordinates.  
// To speed up the pixel lookup, the following function provides a
// matrix which converts output pixel indices to input pixel indices.

vtkMatrix4x4 *vtkImageReslice::GetIndexMatrix()
{
  // first verify that we have to update the matrix
  if (this->IndexMatrix == NULL)
    {
    this->IndexMatrix = vtkMatrix4x4::New();
    }

  int i;
  int isIdentity = 0;
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

  if (this->ResliceAxes)
    {
    transform->SetMatrix(this->GetResliceAxes());
    }
  if (this->ResliceTransform && 
      this->ResliceTransform->IsA("vtkHomogeneousTransform") && 
      this->Optimization)
    {
    transform->PostMultiply();
    transform->Concatenate(((vtkHomogeneousTransform *)
			    this->ResliceTransform)->GetMatrix());
    }
  
  // check to see if we have an identity matrix
  isIdentity = vtkIsIdentityMatrix(transform->GetMatrix());

  // the outMatrix takes OutputData indices to OutputData coordinates,
  // the inMatrix takes InputData coordinates to InputData indices
  for (i = 0; i < 3; i++) 
    {
    if (inSpacing[i] != outSpacing[i] || inOrigin[i] != outOrigin[i])
      {
      isIdentity = 0;
      }
    inMatrix->Element[i][i] = 1.0f/inSpacing[i];
    inMatrix->Element[i][3] = -inOrigin[i]/inSpacing[i];
    outMatrix->Element[i][i] = outSpacing[i];
    outMatrix->Element[i][3] = outOrigin[i];
    };
  this->GetOutput()->GetOrigin(outOrigin); 

  if (!isIdentity)
    {
    transform->PreMultiply();
    transform->Concatenate(outMatrix);
    transform->PostMultiply();
    transform->Concatenate(inMatrix);
    }

  transform->GetMatrix(this->IndexMatrix);
  
  transform->Delete();
  inMatrix->Delete();
  outMatrix->Delete();

  return this->IndexMatrix;
}

//----------------------------------------------------------------------------
// This method is passed a input and output region, and executes the filter
// algorithm to fill the output from the input.
// It just executes a switch statement to call the correct function for
// the regions data types.
void vtkImageReslice::OptimizedThreadedExecute(vtkImageData *inData, 
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

  // change transform matrix so that instead of taking 
  // input coords -> output coords it takes output indices -> input indices
  vtkMatrix4x4 *matrix = this->IndexMatrix;
  
  if (vtkIsPermutationMatrix(matrix))
    {
    switch (inData->GetScalarType())
      {
      case VTK_FLOAT:
	vtkOptimizedPermuteExecute(this, inData, (float *)(inPtr), 
			    outData, (float *)(outPtr),outExt, id,
			    matrix);
	break;
      case VTK_INT:
	vtkOptimizedPermuteExecute(this, inData, (int *)(inPtr), 
			    outData, (int *)(outPtr),outExt, id,
			    matrix);
	break;
      case VTK_SHORT:
	vtkOptimizedPermuteExecute(this, inData, (short *)(inPtr), 
			    outData, (short *)(outPtr),outExt, id,
			    matrix);
	break;
      case VTK_UNSIGNED_SHORT:
	vtkOptimizedPermuteExecute(this, inData, (unsigned short *)(inPtr), 
			    outData, (unsigned short *)(outPtr),outExt,id,
			    matrix);
	break;
      case VTK_UNSIGNED_CHAR:
	vtkOptimizedPermuteExecute(this, inData, (unsigned char *)(inPtr), 
			    outData, (unsigned char *)(outPtr),outExt, id,
			    matrix);
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
	vtkOptimizedExecute(this, inData, (float *)(inPtr), 
			    outData, (float *)(outPtr),outExt, id,
			    matrix);
	break;
      case VTK_INT:
	vtkOptimizedExecute(this, inData, (int *)(inPtr), 
			    outData, (int *)(outPtr),outExt, id,
			    matrix);
	break;
      case VTK_SHORT:
	vtkOptimizedExecute(this, inData, (short *)(inPtr), 
			    outData, (short *)(outPtr),outExt, id,
			    matrix);
	break;
      case VTK_UNSIGNED_SHORT:
	vtkOptimizedExecute(this, inData, (unsigned short *)(inPtr), 
			    outData, (unsigned short *)(outPtr),outExt,id,
			    matrix);
	break;
      case VTK_UNSIGNED_CHAR:
	vtkOptimizedExecute(this, inData, (unsigned char *)(inPtr), 
			    outData, (unsigned char *)(outPtr),outExt, id,
			    matrix);
	break;
      default:
	vtkErrorMacro(<< "Execute: Unknown input ScalarType");
	return;
      }
    }
}






