/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkMarchingContourFilter.cxx
  Language:  C++
  Date:      $Date$
  Version:   $Revision$


Copyright (c) 1993-1998 Ken Martin, Will Schroeder, Bill Lorensen.

    THIS CLASS IS PATENTED UNDER UNITED STATES PATENT NUMBER 4,710,876
    "System and Method for the Display of Surface Structures Contained
    Within the Interior Region of a Solid Body".
    Application of this software for commercial purposes requires 
    a license grant from GE. Contact:
        Jerald Roehling
        GE Licensing
        One Independence Way
        PO Box 2023
        Princeton, N.J. 08540
        phone 609-734-9823
        e-mail:Roehlinj@gerlmo.ge.com
    for more information.

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
#include "vtkMarchingContourFilter.h"
#include "vtkScalars.h"
#include "vtkCell.h"
#include "vtkMergePoints.h"
#include "vtkContourValues.h"
#include "vtkScalarTree.h"

#include "vtkContourFilter.h"
#include "vtkMarchingSquares.h"
#include "vtkMarchingCubes.h"
#include "vtkImageMarchingCubes.h"
#include "vtkObjectFactory.h"



//------------------------------------------------------------------------------
vtkMarchingContourFilter* vtkMarchingContourFilter::New()
{
  // First try to create the object from the vtkObjectFactory
  vtkObject* ret = vtkObjectFactory::CreateInstance("vtkMarchingContourFilter");
  if(ret)
    {
    return (vtkMarchingContourFilter*)ret;
    }
  // If the factory was unable to create the object, then create it here.
  return new vtkMarchingContourFilter;
}




// Construct object with initial range (0,1) and single contour value
// of 0.0.
vtkMarchingContourFilter::vtkMarchingContourFilter()
{
  this->ContourValues = vtkContourValues::New();

  this->ComputeNormals = 1;
  this->ComputeGradients = 0;
  this->ComputeScalars = 1;

  this->Locator = NULL;

  this->UseScalarTree = 0;
  this->ScalarTree = NULL;
}

vtkMarchingContourFilter::~vtkMarchingContourFilter()
{
  this->ContourValues->Delete();
  if ( this->Locator )
    {
    this->Locator->UnRegister(this);
    this->Locator = NULL;
    }
  if ( this->ScalarTree )
    {
    this->ScalarTree->Delete();
    }
}

// Overload standard modified time function. If contour values are modified,
// then this object is modified as well.
unsigned long vtkMarchingContourFilter::GetMTime()
{
  unsigned long mTime=this->vtkDataSetToPolyDataFilter::GetMTime();
  unsigned long time;

  if (this->ContourValues)
    {
    time = this->ContourValues->GetMTime();
    mTime = ( time > mTime ? time : mTime );
    }
  if (this->Locator)
    {
    time = this->Locator->GetMTime();
    mTime = ( time > mTime ? time : mTime );
    }

  return mTime;
}

//
// General contouring filter.  Handles arbitrary input.
//
void vtkMarchingContourFilter::Execute()
{
  vtkScalars *inScalars;
  vtkDataSet *input=this->GetInput();
  vtkPolyData *output=this->GetOutput();
  int numCells;
  vtkPointData *inPd=input->GetPointData(), *outPd=output->GetPointData();
  vtkCellData *inCd=input->GetCellData(), *outCd=output->GetCellData();
  int numContours=this->ContourValues->GetNumberOfContours();
  float *values=this->ContourValues->GetValues();
  
  vtkDebugMacro(<< "Executing marching contour filter");

  numCells = input->GetNumberOfCells();
  inScalars = input->GetPointData()->GetScalars();
  if ( ! inScalars || numCells < 1 )
    {
    vtkErrorMacro(<<"No data to contour");
    return;
    }

  // If structured points, use more efficient algorithms
  if ( (input->GetDataObjectType() == VTK_STRUCTURED_POINTS))
    {
    if (inScalars->GetDataType() != VTK_BIT)
      {
      int dim = input->GetCell(0)->GetCellDimension();
      
      if ( input->GetCell(0)->GetCellDimension() >= 2 ) 
	{
	vtkDebugMacro(<< "Structured Points");
	this->StructuredPointsContour(dim);
	return;
	}
      }
    }
  
  if ( (input->GetDataObjectType() == VTK_IMAGE_DATA)) 
    {
    if (inScalars->GetDataType() != VTK_BIT)
      {
      int dim = input->GetCell(0)->GetCellDimension();
      
      if ( input->GetCell(0)->GetCellDimension() >= 2 ) 
	{
	vtkDebugMacro(<< "Image");
	this->ImageContour(dim);
	return;
	}
      }
    }
  
  vtkDebugMacro(<< "Unoptimized");
  this->DataSetContour();
}

void vtkMarchingContourFilter::StructuredPointsContour(int dim)
{
  vtkPolyData *output;
  vtkPolyData *thisOutput = this->GetOutput();
  int numContours=this->ContourValues->GetNumberOfContours();
  float *values=this->ContourValues->GetValues();

  if ( dim == 2 ) //marching squares
    {
    vtkMarchingSquares *msquares;
    int i;
    
    msquares = vtkMarchingSquares::New();
    msquares->SetInput((vtkStructuredPoints *)this->GetInput());
    msquares->SetDebug(this->Debug);
    msquares->SetNumberOfContours(numContours);
    for (i=0; i < numContours; i++)
      {
      msquares->SetValue(i,values[i]);
      }
         
    msquares->Update();
    output = msquares->GetOutput();
    output->Register(this);
    msquares->Delete();
    }

  else //marching cubes
    {
    vtkMarchingCubes *mcubes;
    int i;
    
    mcubes = vtkMarchingCubes::New();
    mcubes->SetInput((vtkStructuredPoints *)this->GetInput());
    mcubes->SetComputeNormals (this->ComputeNormals);
    mcubes->SetComputeGradients (this->ComputeGradients);
    mcubes->SetComputeScalars (this->ComputeScalars);
    mcubes->SetDebug(this->Debug);
    mcubes->SetNumberOfContours(numContours);
    for (i=0; i < numContours; i++)
      {
      mcubes->SetValue(i,values[i]);
      }

    mcubes->Update();
    output = mcubes->GetOutput();
    output->Register(this);
    mcubes->Delete();
    }
  
  thisOutput->CopyStructure(output);
  thisOutput->GetPointData()->ShallowCopy(output->GetPointData());
  output->UnRegister(this);
}

void vtkMarchingContourFilter::DataSetContour()
{
  vtkPolyData *output = this->GetOutput();
  vtkPolyData *thisOutput = this->GetOutput();
  int numContours=this->ContourValues->GetNumberOfContours();
  float *values=this->ContourValues->GetValues();

  vtkContourFilter *contour = vtkContourFilter::New();
  contour->SetInput((vtkImageData *)this->GetInput());
  contour->SetOutput(output);
  contour->SetComputeNormals (this->ComputeNormals);
  contour->SetComputeGradients (this->ComputeGradients);
  contour->SetComputeScalars (this->ComputeScalars);
  contour->SetDebug(this->Debug);
  contour->SetNumberOfContours(numContours);
  for (int i=0; i < numContours; i++)
    {
    contour->SetValue(i,values[i]);
    }

  contour->Update();
  this->SetOutput(output);
  contour->Delete();
}

void vtkMarchingContourFilter::ImageContour(int dim)
{
  vtkPolyData *output = this->GetOutput();
  vtkPolyData *thisOutput = this->GetOutput();
  int numContours=this->ContourValues->GetNumberOfContours();
  float *values=this->ContourValues->GetValues();

  if ( dim == 2 ) //marching squares
    {
    vtkMarchingSquares *msquares;
    int i;
    
    msquares = vtkMarchingSquares::New();
    msquares->SetInput((vtkImageData *)this->GetInput());
    msquares->SetOutput(output);
    msquares->SetDebug(this->Debug);
    msquares->SetNumberOfContours(numContours);
    for (i=0; i < numContours; i++)
      {
      msquares->SetValue(i,values[i]);
      }
         
    msquares->Update();
    this->SetOutput(output);
    msquares->Delete();
    }

  else //image marching cubes
    {
    vtkImageMarchingCubes *mcubes;
    int i;
    
    mcubes = vtkImageMarchingCubes::New();
    mcubes->SetInput((vtkImageData *)this->GetInput());
    mcubes->SetOutput(output);
    mcubes->SetComputeNormals (this->ComputeNormals);
    mcubes->SetComputeGradients (this->ComputeGradients);
    mcubes->SetComputeScalars (this->ComputeScalars);
    mcubes->SetDebug(this->Debug);
    mcubes->SetNumberOfContours(numContours);
    for (i=0; i < numContours; i++)
      {
      mcubes->SetValue(i,values[i]);
      }

    mcubes->Update();
    this->SetOutput(output);
    mcubes->Delete();
    }
}

// Specify a spatial locator for merging points. By default, 
// an instance of vtkMergePoints is used.
void vtkMarchingContourFilter::SetLocator(vtkPointLocator *locator)
{
  if ( this->Locator == locator ) 
    {
    return;
    }
  if ( this->Locator )
    {
    this->Locator->UnRegister(this);
    this->Locator = NULL;
    }
  if ( locator )
    {
    locator->Register(this);
    }
  this->Locator = locator;
  this->Modified();
}

void vtkMarchingContourFilter::CreateDefaultLocator()
{
  if ( this->Locator == NULL )
    {
    this->Locator = vtkMergePoints::New();
    }
}

void vtkMarchingContourFilter::PrintSelf(ostream& os, vtkIndent indent)
{
  vtkDataSetToPolyDataFilter::PrintSelf(os,indent);

  os << indent << "Compute Gradients: " << (this->ComputeGradients ? "On\n" : "Off\n");
  os << indent << "Compute Normals: " << (this->ComputeNormals ? "On\n" : "Off\n");
  os << indent << "Compute Scalars: " << (this->ComputeScalars ? "On\n" : "Off\n");
  os << indent << "Use Scalar Tree: " << (this->UseScalarTree ? "On\n" : "Off\n");

  this->ContourValues->PrintSelf(os,indent);

  if ( this->Locator )
    {
    os << indent << "Locator: " << this->Locator << "\n";
    }
  else
    {
    os << indent << "Locator: (none)\n";
    }
}

