/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkImageToIsoSurface.cxx
  Language:  C++
  Date:      $Date$
  Version:   $Revision$
  Thanks:    Thanks to C. Charles Law who developed this class.

Copyright (c) 1993-1996 Ken Martin, Will Schroeder, Bill Lorensen.

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
#include "vtkMarchingCubesCases.h"
#include "vtkImageToIsoSurface.h"
#include "vtkImageRegion.h"

//----------------------------------------------------------------------------
vtkImageToIsoSurface::vtkImageToIsoSurface()
{
  for (int i=0; i<VTK_MAX_CONTOURS; i++) this->Values[i] = 0.0;
  this->Input = NULL;
  this->NumberOfContours = 1;
  this->Range[0] = 0.0;
  this->Range[1] = 1.0;
  this->ComputeScalars = 0;
  this->ComputeNormals = 0;
  this->ComputeGradients = 0;

  this->LocatorPointIds = NULL;
  this->InputMemoryLimit = 10000;  // 10 mega Bytes
}


//----------------------------------------------------------------------------
void vtkImageToIsoSurface::PrintSelf(ostream& os, vtkIndent indent)
{
  int i;

  vtkPolySource::PrintSelf(os,indent);

  os << indent << "Number Of Contours : " << this->NumberOfContours << "\n";
  os << indent << "Contour Values: \n";
  for ( i=0; i<this->NumberOfContours; i++)
    {
    os << indent << "  Value " << i << ": " << this->Values[i] << "\n";
    }

  os << indent << "ComputeScalars: " << this->ComputeScalars << "\n";
  os << indent << "ComputeNormals: " << this->ComputeNormals << "\n";
  os << indent << "ComputeGradients: " << this->ComputeGradients << "\n";

  os << indent << "InputMemoryLimit: " << this->InputMemoryLimit <<"K bytes\n";
}


//----------------------------------------------------------------------------
// Description:
// Set a particular contour value at contour number i.
void vtkImageToIsoSurface::SetValue(int i, float value)
{
  i = (i >= VTK_MAX_CONTOURS ? VTK_MAX_CONTOURS-1 : (i < 0 ? 0 : i) );
  if ( this->Values[i] != value )
    {
    this->Modified();
    this->Values[i] = value;
    if ( i >= this->NumberOfContours ) this->NumberOfContours = i + 1;
    if ( value < this->Range[0] ) this->Range[0] = value;
    if ( value > this->Range[1] ) this->Range[1] = value;
    }
}

//----------------------------------------------------------------------------
// Description:
// Generate numContours equally spaced contour values between specified
// range.
void vtkImageToIsoSurface::GenerateValues(int numContours, float range[2])
{
  float val, incr;
  int i;

  numContours = (numContours > VTK_MAX_CONTOURS ? VTK_MAX_CONTOURS : 
                 (numContours > 1 ? numContours : 2) );

  incr = (range[1] - range[0]) / (numContours-1);
  for (i=0, val=range[0]; i < numContours; i++, val+=incr)
    {
    this->SetValue(i,val);
    }
}

//----------------------------------------------------------------------------
// Description:
// Generate numContours equally spaced contour values between specified
// range.
void vtkImageToIsoSurface::GenerateValues(int numContours, float r1, float r2)
{
  float rng[2];

  rng[0] = r1;
  rng[1] = r2;
  this->GenerateValues(numContours,rng);
}

//----------------------------------------------------------------------------
void vtkImageToIsoSurface::Execute()
{
  vtkImageRegion *inRegion;
  vtkPolyData *output = this->GetOutput();
  int extent[6], estimatedSize;
  int temp, zMin, zMax, chunkMin, chunkMax;
  int minSlicesPerChunk, chunkOverlap;
  
  if ( ! this->Input)
    {
    vtkErrorMacro(<< "No Input");
    return;
    }
  
  // Gradients must be computed (but not saved) if Compute normals is on.
  this->NeedGradients = this->ComputeGradients || this->ComputeNormals;
  
  // Determine the number of slices per request from input memory limit.
  if (this->NeedGradients)
    {
    minSlicesPerChunk = 4;
    chunkOverlap = 3;
    }
  else
    {
    minSlicesPerChunk = 2;
    chunkOverlap = 1;
    }
  inRegion = vtkImageRegion::New();
  this->Input->UpdateImageInformation(inRegion);
  // Each data type requires a different amount of memory.
  switch (this->Input->GetScalarType())
    {
    case VTK_FLOAT:
      temp = sizeof(float);
      break;
    case VTK_INT:
      temp = sizeof(int);
      break;
    case VTK_SHORT:
      temp = sizeof(short);
      break;
    case VTK_UNSIGNED_SHORT:
      temp = sizeof(unsigned short);
      break;
    case VTK_UNSIGNED_CHAR:
      temp = sizeof(unsigned char);
      break;
    default:
      vtkErrorMacro(<< "Could not determine input scalar type.");
      return;
    }
  inRegion->GetImageExtent(3, extent);
  // multiply by the area of each slice
  temp *= extent[1] - extent[0] + 1;
  temp *= extent[3] - extent[2] + 1;
  temp = temp;
  // temp holds memory per image. (+1 to avoid dividing by zero)
  this->NumberOfSlicesPerChunk = this->InputMemoryLimit * 1000 / (temp + 1);
  if (this->NumberOfSlicesPerChunk < minSlicesPerChunk)
    {
    vtkWarningMacro("Execute: Need " << (temp/1000) << " KB to load " 
		    << minSlicesPerChunk << " minimum.\n");
    this->NumberOfSlicesPerChunk = minSlicesPerChunk;
    }
  vtkDebugMacro("Execute: NumberOfSlicesPerChunk = " 
		<< this->NumberOfSlicesPerChunk);
  this->NumberOfSlicesPerChunk -= chunkOverlap;
  
  // Create the points, scalars, normals and Cell arrays for the output.
  // estimate the number of points from the volume dimensions
  estimatedSize = (int) pow ((double) ((extent[1]-extent[0]+1) * 
				       (extent[3]-extent[2]+1) * 
				       (extent[5]-extent[4]+1)), 0.75);
  estimatedSize = estimatedSize / 1024 * 1024; //multiple of 1024
  if (estimatedSize < 1024) 
    {
    estimatedSize = 1024;
    }
  vtkDebugMacro(<< "Estimated number of points/triangles: " << estimatedSize);
  this->Points = new vtkFloatPoints(estimatedSize,estimatedSize/2);
  this->Triangles = new vtkCellArray(estimatedSize,estimatedSize/2);
  if (this->ComputeScalars)
    {
    this->Scalars = new vtkFloatScalars(estimatedSize,estimatedSize/2);
    }
  if (this->ComputeNormals)
    {
    this->Normals = new vtkFloatNormals(estimatedSize,estimatedSize/2);
    }
  if (this->ComputeGradients)
    {
    this->Gradients = new vtkFloatVectors(estimatedSize,estimatedSize/2);
    }

  // Initialize the internal point locator (edge table for oen image of cubes).
  this->InitializeLocator(extent[0], extent[1], extent[2], extent[3]);
  
  // Loop through the chunks running marching cubes on each one
  zMin = extent[4];
  zMax = extent[5];
  for(chunkMin = zMin; chunkMin < zMax; chunkMin = chunkMax)
    {
    // Get the chunk from the input
    chunkMax = chunkMin + this->NumberOfSlicesPerChunk;
    if (chunkMax > zMax)
      {
      chunkMax = zMax;
      }
    extent[4] = chunkMin;
    extent[5] = chunkMax;
    // Expand if computing gradients with central differences
    if (this->NeedGradients)
      {
      --extent[4];
      ++extent[5];
      }
    // Don't go over boundary of data.
    if (extent[4] < zMin)
      {
      extent[4] = zMin;
      }
    if (extent[5] > zMax)
      {
      extent[5] = zMax;
      }
    inRegion->SetExtent(3, extent);
    // Get the chunk from the input
    this->Input->UpdateRegion(inRegion);
    
    this->March(inRegion, chunkMin, chunkMax);
    }
  
  // Put results in our output
  vtkDebugMacro(<<"Created: " 
               << this->Points->GetNumberOfPoints() << " points, " 
               << this->Triangles->GetNumberOfCells() << " triangles");
  output->SetPoints(this->Points);
  this->Points->Delete();
  this->Points = NULL;
  output->SetPolys(this->Triangles);
  this->Triangles->Delete();
  this->Triangles = NULL;
  if (this->ComputeScalars)
    {
    output->GetPointData()->SetScalars(this->Scalars);
    this->Scalars->Delete();
    this->Scalars = NULL;
    }
  if (this->ComputeNormals)
    {
    output->GetPointData()->SetNormals(this->Normals);
    this->Normals->Delete();
    this->Normals = NULL;
    }
  
  // Recover extra space.
  output->Squeeze();
  
  // Clean up temporary images.
  inRegion->Delete();
  // release the locators memory
  this->DeleteLocator();
}





//----------------------------------------------------------------------------
// This method selectively applies marching cubes (iso surface = 0.0)
// to the derivative (second derivative because vector was first).
// the cube is ignored if the magnitude values are below the 
// MagnitudeThreshold.
void vtkImageToIsoSurface::March(vtkImageRegion *inRegion, 
				 int chunkMin, int chunkMax)
{
  int idx0, idx1, idx2;
  int min0, max0, min1, max1, min2, max2;
  int inc0, inc1, inc2;
  float *ptr0, *ptr1, *ptr2;
  
  // Get information to loop through images.
  inRegion->GetExtent(min0, max0, min1, max1, min2, max2);
  ptr2 = (float *)(inRegion->GetScalarPointer(min0, min1, chunkMin));
  inRegion->GetIncrements(inc0, inc1, inc2);

  // Loop over all the cubes
  for (idx2 = chunkMin; idx2 < chunkMax; ++idx2)
    {
    ptr1 = ptr2;
    for (idx1 = min1; idx1 < max1; ++idx1)
      {
      ptr0 = ptr1;
      for (idx0 = min0; idx0 < max0; ++idx0)
	{
	// put magnitudes into the cube structure.
	this->HandleCube(idx0, idx1, idx2, inRegion, ptr0);

	ptr0 += inc0;
	}
      ptr1 += inc1;
      }
    ptr2 += inc2;
    this->IncrementLocatorZ();
    }
}





//----------------------------------------------------------------------------
// This method runs marching cubes on one cube.
void vtkImageToIsoSurface::HandleCube(int cellX, int cellY, int cellZ,
				      vtkImageRegion *inRegion,
				      float *ptr)
{
  int inc0, inc1, inc2;
  int valueIdx;
  float value;
  int cubeIndex, ii, pointIds[3];
  TRIANGLE_CASES *triCase;
  EDGE_LIST  *edge;

  inRegion->GetIncrements(inc0, inc1, inc2);
  for (valueIdx = 0; valueIdx < this->NumberOfContours; ++valueIdx)
    {
    value = this->Values[valueIdx];
    // compute the case index
    cubeIndex = 0;
    if (ptr[0] > value) cubeIndex += 1;
    if (ptr[inc0] > value) cubeIndex += 2;
    if (ptr[inc0 + inc1] > value) cubeIndex += 4;
    if (ptr[inc1] > value) cubeIndex += 8;
    if (ptr[inc2] > value) cubeIndex += 16;
    if (ptr[inc0 + inc2] > value) cubeIndex += 32;
    if (ptr[inc0 + inc1 + inc2] > value) cubeIndex += 64;
    if (ptr[inc1 + inc2] > value) cubeIndex += 128;
    // Make sure we have trianlges
    if (cubeIndex != 0 && cubeIndex != 255)
      {
      // Get edges.
      triCase = triCases + cubeIndex;
      edge = triCase->edges; 
      // loop over triangles  
      while(*edge > -1)
	{
	for (ii=0; ii<3; ++ii, ++edge) //insert triangle
	  {
	  // Get the index of the point
	  pointIds[ii] = this->GetLocatorPoint(cellX, cellY, *edge);
	  // If the point has not been created yet
	  if (pointIds[ii] == -1)
	    {
	    float *aspectRatio = inRegion->GetAspectRatio();
	    float *origin = inRegion->GetOrigin();
	    int *extent = inRegion->GetImageExtent();
	    
	    pointIds[ii] = this->MakeNewPoint(cellX, cellY, cellZ,
					      inc0, inc1, inc2,
					      ptr, *edge, extent,
					      aspectRatio, origin, value);
	    this->AddLocatorPoint(cellX, cellY, *edge, pointIds[ii]);
	    }
	  }
	this->Triangles->InsertNextCell(3,pointIds);
	}//for each triangle
      }
    }
}





//----------------------------------------------------------------------------
// This method interpolates verticies to make a new point.
int vtkImageToIsoSurface::MakeNewPoint(int idx0, int idx1, int idx2,
				       int inc0, int inc1, int inc2,
				       float *ptr, int edge, int *imageExtent,
				       float *aspectRatio, float *origin,
				       float value)
{
  int edgeAxis;
  float *ptrB;
  float temp, pt[3];

  // decode the edge into starting point and axis direction
  switch (edge)
    {
    case 0:  // 0,1
      ptrB = ptr + inc0;
      edgeAxis = 0;
      break;
    case 1:  // 1,2
      ++idx0;
      ptr += inc0;
      ptrB = ptr + inc1;
      edgeAxis = 1;
      break;
    case 2:  // 3,2
      ++idx1;
      ptr += inc1;
      ptrB = ptr + inc0;
      edgeAxis = 0;
      break;
    case 3:  // 0,3
      ptrB = ptr + inc1;
      edgeAxis = 1;
      break;
    case 4:  // 4,5
      ++idx2;
      ptr += inc2;
      ptrB = ptr + inc0;
      edgeAxis = 0;
      break;
    case 5:  // 5,6
      ++idx0; ++idx2;
      ptr += inc0 + inc2;
      ptrB = ptr + inc1;
      edgeAxis = 1;
      break;
    case 6:  // 7,6
      ++idx1; ++idx2;
      ptr += inc1 + inc2;
      ptrB = ptr + inc0;
      edgeAxis = 0;
      break;
    case 7: // 4,7
      ++idx2;
      ptr += inc2;
      ptrB = ptr + inc1;
      edgeAxis = 1;
      break;
    case 8: // 0,4
      ptrB = ptr + inc2;
      edgeAxis = 2;
      break;
    case 9: // 1,5
      ++idx0;
      ptr += inc0;
      ptrB = ptr + inc2;
      edgeAxis = 2;
      break;
    case 10: // 3,7
      ++idx1;
      ptr += inc1;
      ptrB = ptr + inc2;
      edgeAxis = 2;
      break;
    case 11: // 2,6
      ++idx0; ++idx1;
      ptr += inc0 + inc1;
      ptrB = ptr + inc2;
      edgeAxis = 2;
      break;
    }

  // interpolation factor
  temp = (value - *ptr) / (*ptrB - *ptr);
  
  // interpolate the point position
  switch (edgeAxis)
    {
    case 0:
      pt[0] = origin[0] + aspectRatio[0] * ((float)idx0 + temp);
      pt[1] = origin[1] + aspectRatio[1] * ((float)idx1);
      pt[2] = origin[2] + aspectRatio[2] * ((float)idx2);
      break;
    case 1:
      pt[0] = origin[0] + aspectRatio[0] * ((float)idx0);
      pt[1] = origin[1] + aspectRatio[1] * ((float)idx1 + temp);
      pt[2] = origin[2] + aspectRatio[2] * ((float)idx2);
      break;
    case 2:
      pt[0] = origin[0] + aspectRatio[0] * ((float)idx0);
      pt[1] = origin[1] + aspectRatio[1] * ((float)idx1);
      pt[2] = origin[2] + aspectRatio[2] * ((float)idx2 + temp);
      break;
    }
  
  // Save the scale if we are generating scalars
  if (this->ComputeScalars)
    {
    this->Scalars->InsertNextScalar(value);
    }
  
  // Interpolate to find normal from vectors.
  if (this->NeedGradients)
    {
    short b0, b1, b2;
    float g[3], gB[3];
    // Find boundary conditions and compute gradient (first point)
    b0 = (idx0 == imageExtent[1]);
    if (idx0 == imageExtent[0]) b0 = -1;
    b1 = (idx1 == imageExtent[3]);
    if (idx1 == imageExtent[2]) b1 = -1;
    b2 = (idx2 == imageExtent[5]);
    if (idx2 == imageExtent[4]) b2 = -1;
    this->ComputePointGradient(g, ptr, inc0, inc1, inc2, b0, b1, b2);
    // Find boundary conditions and compute gradient (second point)
    switch (edgeAxis)
      {
      case 0:  
	++idx0;
	b0 = (idx0 == imageExtent[1]);
	break;
      case 1:  
	++idx1;
	b1 = (idx1 == imageExtent[3]);
	break;
      case 2:  
	++idx2;
	b2 = (idx2 == imageExtent[5]);
	break;
      }
    this->ComputePointGradient(gB, ptrB, inc0, inc1, inc2, b0, b1, b2);
    // Interpolate Gradient
    g[0] = g[0] + temp * (gB[0] - g[0]);
    g[1] = g[1] + temp * (gB[1] - g[1]);
    g[2] = g[2] + temp * (gB[2] - g[2]);
    if (this->ComputeGradients)
      {
      this->Gradients->InsertNextVector(g);
      }
    if (this->ComputeNormals)
      {
      temp = -1.0 / sqrt(g[0]*g[0] + g[1]*g[1] + g[2]*g[2]);
      g[0] *= temp;
      g[1] *= temp;
      g[2] *= temp;
      this->Normals->InsertNextNormal(g);
      }
    }
  
  return this->Points->InsertNextPoint(pt);
}


//----------------------------------------------------------------------------
// This method uses central differences to compute the gradient
// of a point. Note: This method assumes that max > min for all 3 axes!
// It does not consider aspect ratio.
void vtkImageToIsoSurface::ComputePointGradient(float g[3], float *ptr,
						int inc0, int inc1, int inc2,
						int b0, int b1, int b2)
{
  if (b0 < 0)
    {
    g[0] = ptr[inc0] - *ptr;
    }
  else if (b0 > 0)
    {
    g[0] = *ptr - ptr[-inc0];
    }
  else
    {
    g[0] = ptr[inc0] - ptr[-inc0];
    }

  if (b1 < 0)
    {
    g[1] = ptr[inc1] - *ptr;
    }
  else if (b1 > 0)
    {
    g[1] = *ptr - ptr[-inc1];
    }
  else
    {
    g[1] = ptr[inc1] - ptr[-inc1];
    }

  if (b2 < 0)
    {
    g[2] = ptr[inc2] - *ptr;
    }
  else if (b2 > 0)
    {
    g[2] = *ptr - ptr[-inc2];
    }
  else
    {
    g[2] = ptr[inc2] - ptr[-inc2];
    }
}


//============================================================================
// These method act as the point locator so verticies will be shared.
// One 2d array of cubes is stored. (z dimension is ignored).
// Points are indexed by their cube and edge.
// Shared edges are only represented once.  Cubes are responsible for
// edges on their min faces.  Their is an extra row and column of cubes
// to store the max edges of the last row/column of cubes,


//----------------------------------------------------------------------------
// This method allocates and initializes the point array.
// One 2d array of cubes is stored. (z dimension is ignored).
void vtkImageToIsoSurface::InitializeLocator(int min0, int max0, 
						 int min1, int max1)
{
  int idx;
  int size;

  // Free old memory
  if (this->LocatorPointIds)
    {
    delete this->LocatorPointIds;
    }
  // Extra row and column
  this->LocatorDimX = (max0 - min0 + 2);
  this->LocatorDimY = (max1 - min1 + 2);
  this->LocatorMinX = min0;
  this->LocatorMinY = min1;
  // 5 non shared edges.
  size = (this->LocatorDimX)*(this->LocatorDimY)*5;
  this->LocatorPointIds = new int[size];
  // Initialize the array
  for (idx = 0; idx < size; ++idx)
    {
    this->LocatorPointIds[idx] = -1;
    }
}

//----------------------------------------------------------------------------
// This method frees the locators memory.
void vtkImageToIsoSurface::DeleteLocator()
{
  // Free old memory
  if (this->LocatorPointIds)
    {
    delete this->LocatorPointIds;
    }
}

//----------------------------------------------------------------------------
// This method moves the Z index of the locator up one slice.
void vtkImageToIsoSurface::IncrementLocatorZ()
{
  int x, y;
  int *ptr;

  ptr = this->LocatorPointIds;
  for (y = 0; y < this->LocatorDimY; ++y)
    {
    for (x = 0; x < this->LocatorDimX; ++x)
      {
      ptr[0] = ptr[4];
      ptr[3] = ptr[1];
      ptr[1] = ptr[2] = ptr[4] = -1;
      ptr += 5;
      }
    }
}

//----------------------------------------------------------------------------
// This method adds a point to the array.  Cube is the X/Y cube,
// segment is the index of the segment (same as marching cubes).(XYZ)
// (0,0,0)->(1,0,0): 0,  (1,0,0)->(1,1,0): 1,
// (1,1,0)->(0,1,0): 2,  (0,1,0)->(0,0,0): 3,
// (0,0,1)->(1,0,1): 4,  (1,0,1)->(1,1,1): 5,
// (1,1,1)->(0,1,1): 6,  (0,1,1)->(0,0,1): 7,
// (0,0,0)->(0,0,1): 8,  (1,0,0)->(1,0,1): 9,
// (0,1,0)->(0,1,1): 10, (1,1,0)->(1,1,1): 11.
// Shared edges are computed internaly. (no error checking)
void vtkImageToIsoSurface::AddLocatorPoint(int cellX, int cellY, int edge,
					       int ptId)
{
  int *ptr;
  
  // Get the correct position in the array.
  ptr = this->GetLocatorPointer(cellX, cellY, edge);
  *ptr = ptId;
}

//----------------------------------------------------------------------------
// This method gets a point from the locator.
int vtkImageToIsoSurface::GetLocatorPoint(int cellX, int cellY, int edge)
{
  int *ptr;
  
  // Get the correct position in the array.
  ptr = this->GetLocatorPointer(cellX, cellY, edge);
  return *ptr;
}

//----------------------------------------------------------------------------
// This method returns a pointer to an ID from a cube and an edge.
int *vtkImageToIsoSurface::GetLocatorPointer(int cellX,int cellY,int edge)
{
  // Remove redundant edges (shared by more than one cube).
  // Take care of shared edges
  switch (edge)
    {
    case 9:  ++cellX;          edge = 8; break;
    case 10: ++cellY;          edge = 8; break;
    case 11: ++cellX; ++cellY; edge = 8; break;
    case 5:  ++cellX;          edge = 7; break;
    case 6:  ++cellY;          edge = 4; break;
    case 1:  ++cellX;          edge = 3; break;
    case 2:  ++cellY;          edge = 0; break;
    }
  
  // relative to min and max.
  cellX -= this->LocatorMinX;
  cellY -= this->LocatorMinY;

  // compute new indexes for edges (0 to 4) 
  // must be compatable with LocatorIncrementZ.
  if (edge == 7)
    {
    edge = 1;
    }
  if (edge == 8)
    {
    edge = 2;
    }

  // return correct pointer
  return this->LocatorPointIds + edge 
    + (cellX + cellY * (this->LocatorDimX)) * 5;
}

  
  





