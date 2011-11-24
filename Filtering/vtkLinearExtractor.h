/*=========================================================================

Program:   Visualization Toolkit
Module:    $RCSfile: vtkCellDistanceFilter,v $

Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
All rights reserved.
See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkLinearExtractor - select cells intersecting a line
//
// .SECTION Description
// This filter takes a vtkCompositeDataSet as input and a line segment as parameter. 
// It outputs a vtkSelection identifying all the cells intersecting the given line segment.
//
// .SECTION Thanks
// This file has been initially developed in the frame of CEA's Love visualization software development <br>
// CEA/DIF - Commissariat a l'Energie Atomique, Centre DAM Ile-De-France <br>
// BP12, F-91297 Arpajon, France. <br>
// This class was implemented by Thierry Carrard, Charles Pignerol, and Philippe Pebay, Kitware, 2011.

#ifndef VTK_LINEAR_EXTRACTOR_H
#define VTK_LINEAR_EXTRACTOR_H

#include <vtkSelectionAlgorithm.h>

class vtkAlgorithmOutput;
class vtkDataSet;
class vtkDoubleArray;
class vtkIdTypeArray;
class vtkPoints;

class VTK_FILTERING_EXPORT vtkLinearExtractor: public vtkSelectionAlgorithm
{
 public:
  vtkTypeMacro(vtkLinearExtractor,vtkSelectionAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  static vtkLinearExtractor *New();

  // Description:
  // Set/Get starting point of intersecting segment
  vtkSetVector3Macro(StartPoint,double);
  vtkGetVectorMacro(StartPoint,double,3);

  // Description:
  // Set/Get end point of intersecting segment
  vtkSetVector3Macro(EndPoint,double);
  vtkGetVectorMacro(EndPoint,double,3);

  // Description:
  // Set/Get the list of points defining the intersecting broken line
  virtual void SetPoints(vtkPoints*);
  vtkGetObjectMacro(Points,vtkPoints);

  // Description:
  // Set/Get tolerance to be used by intersection algorithm
  vtkSetMacro(Tolerance,double);
  vtkGetMacro(Tolerance,double);

  // Description:
  // Set/Get whether lines vertice are included in selection
  vtkSetMacro( IncludeVertices, bool );
  vtkGetMacro( IncludeVertices, bool );

  // Description:
  // Set/Get relative tolerance for vertex elimination
  vtkSetClampMacro(VertexEliminationTolerance,double,0.,.1 );
  vtkGetMacro(VertexEliminationTolerance,double);

 protected:
  vtkLinearExtractor();
  virtual ~vtkLinearExtractor();

  virtual int FillInputPortInformation(int port, vtkInformation *info);

  virtual int RequestData(vtkInformation *request,
                          vtkInformationVector **inputVector,
                          vtkInformationVector *outputVector);

  void RequestDataInternal(vtkDataSet* input, vtkIdTypeArray* outIndices);

 private:
  vtkLinearExtractor(const vtkLinearExtractor&);  // Not implemented
  void operator =(const vtkLinearExtractor&); // Not implemented

  // Description:
  // Start and end point of the intersecting line segment
  // NB: These are used if and only if Points is null.
  double StartPoint[3];
  double EndPoint[3];

  // Description:
  // The list of points defining the intersecting broken line
  // NB: The Start/EndPoint definition of a single line segment is used by default
  vtkPoints* Points;

  // Description:
  // Tolerance to be used by intersection algorithm
  double Tolerance;

  // Description:
  // Decide whether lines vertice are included in selection
  // Default: true
  bool IncludeVertices;

  // Description:
  // Relative tolerance for vertex elimination
  // Default: 1e-6
  double VertexEliminationTolerance;
};


#endif	// VTK_LINEAR_EXTRACTOR_H
