/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkMarchingContourFilter.h
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

  Copyright (c) 1993-2002 Ken Martin, Will Schroeder, Bill Lorensen 
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notice for more information.

     THIS CLASS IS PATENTED UNDER UNITED STATES PATENT NUMBER 4,710,876
     "System and Method for the Display of Surface Structures Contained
     Within the Interior Region of a Solid Body".
     Application of this software for commercial purposes requires 
     a license grant from GE. Contact:

         Carl B. Horton
         Sr. Counsel, Intellectual Property
         3000 N. Grandview Blvd., W-710
         Waukesha, WI  53188
         Phone:  (262) 513-4022
         E-Mail: Carl.Horton@med.ge.com

     for more information.

=========================================================================*/
// .NAME vtkMarchingContourFilter - generate isosurfaces/isolines from scalar values
// .SECTION Description
// vtkMarchingContourFilter is a filter that takes as input any dataset and 
// generates on output isosurfaces and/or isolines. The exact form 
// of the output depends upon the dimensionality of the input data. 
// Data consisting of 3D cells will generate isosurfaces, data 
// consisting of 2D cells will generate isolines, and data with 1D 
// or 0D cells will generate isopoints. Combinations of output type 
// are possible if the input dimension is mixed.
//
// This filter will identify special dataset types (e.g., structured
// points) and use the appropriate specialized filter to process the
// data. For examples, if the input dataset type is a volume, this
// filter will create an internal vtkMarchingCubes instance and use
// it. This gives much better performance.
// 
// To use this filter you must specify one or more contour values.
// You can either use the method SetValue() to specify each contour
// value, or use GenerateValues() to generate a series of evenly
// spaced contours. It is also possible to accelerate the operation of
// this filter (at the cost of extra memory) by using a
// vtkScalarTree. A scalar tree is used to quickly locate cells that
// contain a contour surface. This is especially effective if multiple
// contours are being extracted. If you want to use a scalar tree,
// invoke the method UseScalarTreeOn().

// .SECTION Caveats
// For unstructured data or structured grids, normals and gradients
// are not computed.  This calculation will be implemented in the
// future. In the mean time, use vtkPolyDataNormals to compute the surface
// normals.

// .SECTION See Also
// vtkMarchingCubes vtkSliceCubes vtkDividingCubes vtkMarchingSquares
// vtkImageMarchingCubes

#ifndef __vtkMarchingContourFilter_h
#define __vtkMarchingContourFilter_h

#include "vtkDataSetToPolyDataFilter.h"
#include "vtkContourValues.h"

class vtkScalarTree;

class VTK_PATENTED_EXPORT vtkMarchingContourFilter : public vtkDataSetToPolyDataFilter
{
public:
  vtkTypeRevisionMacro(vtkMarchingContourFilter,vtkDataSetToPolyDataFilter);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Construct object with initial range (0,1) and single contour value
  // of 0.0.
  static vtkMarchingContourFilter *New();

  // Description:
  // Methods to set / get contour values.
  void SetValue(int i, float value);
  float GetValue(int i);
  float *GetValues();
  void GetValues(float *contourValues);
  void SetNumberOfContours(int number);
  int GetNumberOfContours();
  void GenerateValues(int numContours, float range[2]);
  void GenerateValues(int numContours, float rangeStart, float rangeEnd);

  // Description:
  // Modified GetMTime Because we delegate to vtkContourValues
  unsigned long GetMTime();

  // Description:
  // Set/Get the computation of normals. Normal computation is fairly
  // expensive in both time and storage. If the output data will be
  // processed by filters that modify topology or geometry, it may be
  // wise to turn Normals and Gradients off.
  vtkSetMacro(ComputeNormals,int);
  vtkGetMacro(ComputeNormals,int);
  vtkBooleanMacro(ComputeNormals,int);

  // Description:
  // Set/Get the computation of gradients. Gradient computation is
  // fairly expensive in both time and storage. Note that if
  // ComputeNormals is on, gradients will have to be calculated, but
  // will not be stored in the output dataset.  If the output data
  // will be processed by filters that modify topology or geometry, it
  // may be wise to turn Normals and Gradients off.
  vtkSetMacro(ComputeGradients,int);
  vtkGetMacro(ComputeGradients,int);
  vtkBooleanMacro(ComputeGradients,int);

  // Description:
  // Set/Get the computation of scalars.
  vtkSetMacro(ComputeScalars,int);
  vtkGetMacro(ComputeScalars,int);
  vtkBooleanMacro(ComputeScalars,int);

  // Description:
  // Enable the use of a scalar tree to accelerate contour extraction.
  vtkSetMacro(UseScalarTree,int);
  vtkGetMacro(UseScalarTree,int);
  vtkBooleanMacro(UseScalarTree,int);

  // Description:
  // Set / get a spatial locator for merging points. By default, 
  // an instance of vtkMergePoints is used.
  void SetLocator(vtkPointLocator *locator);
  vtkGetObjectMacro(Locator,vtkPointLocator);

  // Description:
  // Create default locator. Used to create one when none is
  // specified. The locator is used to merge coincident points.
  void CreateDefaultLocator();

protected:
  vtkMarchingContourFilter();
  ~vtkMarchingContourFilter();

  void Execute();

  vtkContourValues *ContourValues;
  int ComputeNormals;
  int ComputeGradients;
  int ComputeScalars;
  vtkPointLocator *Locator;
  int UseScalarTree;
  vtkScalarTree *ScalarTree;
  
  //special contouring for structured points
  void StructuredPointsContour(int dim); 
  //special contouring for image data
  void ImageContour(int dim);
  //default if not structured data
  void DataSetContour();
private:
  vtkMarchingContourFilter(const vtkMarchingContourFilter&);  // Not implemented.
  void operator=(const vtkMarchingContourFilter&);  // Not implemented.
};

// Description:
// Set a particular contour value at contour number i. The index i ranges 
// between 0<=i<NumberOfContours.
inline void vtkMarchingContourFilter::SetValue(int i, float value)
{
  this->ContourValues->SetValue(i,value);
}

// Description:
// Get the ith contour value.
inline float vtkMarchingContourFilter::GetValue(int i)
{
  return this->ContourValues->GetValue(i);
}

// Description:
// Get a pointer to an array of contour values. There will be
// GetNumberOfContours() values in the list.
inline float *vtkMarchingContourFilter::GetValues()
{
  return this->ContourValues->GetValues();
}

// Description:
// Fill a supplied list with contour values. There will be
// GetNumberOfContours() values in the list. Make sure you allocate
// enough memory to hold the list.
inline void vtkMarchingContourFilter::GetValues(float *contourValues)
{
  this->ContourValues->GetValues(contourValues);
}

// Description:
// Set the number of contours to place into the list. You only really
// need to use this method to reduce list size. The method SetValue()
// will automatically increase list size as needed.
inline void vtkMarchingContourFilter::SetNumberOfContours(int number)
{
  this->ContourValues->SetNumberOfContours(number);
}

// Description:
// Get the number of contours in the list of contour values.
inline int vtkMarchingContourFilter::GetNumberOfContours()
{
  return this->ContourValues->GetNumberOfContours();
}

// Description:
// Generate numContours equally spaced contour values between specified
// range. Contour values will include min/max range values.
inline void vtkMarchingContourFilter::GenerateValues(int numContours,
                                                     float range[2])
{
  this->ContourValues->GenerateValues(numContours, range);
}

// Description:
// Generate numContours equally spaced contour values between specified
// range. Contour values will include min/max range values.
inline void vtkMarchingContourFilter::GenerateValues(int numContours,
                                                     float rangeStart,
                                                     float rangeEnd)
{
  this->ContourValues->GenerateValues(numContours, rangeStart, rangeEnd);
}

#endif
