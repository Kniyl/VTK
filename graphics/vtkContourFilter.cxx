/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkContourFilter.cxx
  Language:  C++
  Date:      $Date$
  Version:   $Revision$


Copyright (c) 1993-1998 Ken Martin, Will Schroeder, Bill Lorensen.

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
#include "vtkContourFilter.h"
#include "vtkScalars.h"
#include "vtkCell.h"
#include "vtkMergePoints.h"
#include "vtkContourValues.h"
#include "vtkScalarTree.h"

#ifdef VTK_USE_PATENTED
#include "vtkMarchingSquares.h"
#include "vtkMarchingCubes.h"
#include "vtkImageMarchingCubes.h"
#endif

// Construct object with initial range (0,1) and single contour value
// of 0.0.
vtkContourFilter::vtkContourFilter()
{
  this->ContourValues = vtkContourValues::New();

  this->ComputeNormals = 1;
  this->ComputeGradients = 0;
  this->ComputeScalars = 1;

  this->Locator = NULL;

  this->UseScalarTree = 0;
  this->ScalarTree = NULL;
}

vtkContourFilter::~vtkContourFilter()
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
unsigned long vtkContourFilter::GetMTime()
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
void vtkContourFilter::Execute()
{
  int cellId, i;
  vtkIdList *cellPts;
  vtkScalars *inScalars;
  vtkCell *cell;
  float range[2];
  vtkCellArray *newVerts, *newLines, *newPolys;
  vtkPoints *newPts;
  vtkDataSet *input=this->GetInput();
  vtkPolyData *output=this->GetOutput();
  int numCells, estimatedSize;
  vtkPointData *inPd=input->GetPointData(), *outPd=output->GetPointData();
  vtkCellData *inCd=input->GetCellData(), *outCd=output->GetCellData();
  int numContours=this->ContourValues->GetNumberOfContours();
  float *values=this->ContourValues->GetValues();
  vtkScalars *cellScalars;
  
  vtkDebugMacro(<< "Executing contour filter");

  numCells = input->GetNumberOfCells();
  inScalars = input->GetPointData()->GetScalars();
  if ( ! inScalars || numCells < 1 )
    {
    vtkErrorMacro(<<"No data to contour");
    return;
    }

  // If structured points, use more efficient algorithms
#ifdef VTK_USE_PATENTED
  if ( input->GetDataSetType() == VTK_STRUCTURED_POINTS &&
       inScalars->GetDataType() != VTK_BIT)
    {
    int dim = input->GetCell(0)->GetCellDimension();

    if ( input->GetCell(0)->GetCellDimension() >= 2 ) 
      {
      this->StructuredPointsContour(dim);
      return;
      }
    }
#endif

  inScalars->GetRange(range);
//
// Create objects to hold output of contour operation. First estimate allocation size.
//
  estimatedSize = (int) pow ((double) numCells, .75);
  estimatedSize *= numContours;
  estimatedSize = estimatedSize / 1024 * 1024; //multiple of 1024
  if (estimatedSize < 1024)
    {
    estimatedSize = 1024;
    }

  newPts = vtkPoints::New();
  newPts->Allocate(estimatedSize,estimatedSize);
  newVerts = vtkCellArray::New();
  newVerts->Allocate(estimatedSize,estimatedSize);
  newLines = vtkCellArray::New();
  newLines->Allocate(estimatedSize,estimatedSize);
  newPolys = vtkCellArray::New();
  newPolys->Allocate(estimatedSize,estimatedSize);
  cellScalars = vtkScalars::New();
  cellScalars->Allocate(VTK_CELL_SIZE);
  
 
  // locator used to merge potentially duplicate points
  if ( this->Locator == NULL )
    {
    this->CreateDefaultLocator();
    }
  this->Locator->InitPointInsertion (newPts, input->GetBounds(),estimatedSize);

  // interpolate data along edge
  // if we did not ask for scalars to be computed, don't copy them
  if (!this->ComputeScalars)
    {
    outPd->CopyScalarsOff();
    }
  outPd->InterpolateAllocate(inPd,estimatedSize,estimatedSize);
  outCd->CopyAllocate(inCd,estimatedSize,estimatedSize);

  // If enabled, build a scalar tree to accelerate search
  //
  if ( !this->UseScalarTree )
    {
    for (cellId=0; cellId < numCells; cellId++)
      {
      cell = input->GetCell(cellId);
      cellPts = cell->GetPointIds();
      inScalars->GetScalars(cellPts,cellScalars);

      for (i=0; i < numContours; i++)
        {
        cell->Contour(values[i], cellScalars, this->Locator, 
                      newVerts, newLines, newPolys, inPd, outPd,
		      inCd, cellId, outCd);

        } // for all contour values
      } // for all cells
    } //if using scalar tree
  else
    {
    if ( this->ScalarTree == NULL )
      {
      this->ScalarTree = new vtkScalarTree;
      }
    this->ScalarTree->SetDataSet(input);
    //
    // Loop over all contour values.  Then for each contour value, 
    // loop over all cells.
    //
    for (i=0; i < numContours; i++)
      {
      for ( this->ScalarTree->InitTraversal(values[i]); 
      (cell=this->ScalarTree->GetNextCell(cellId,cellPts,cellScalars)) != NULL; )
        {
        cell->Contour(values[i], cellScalars, this->Locator, 
                      newVerts, newLines, newPolys, inPd, outPd,
		      inCd, cellId, outCd);

        } //for all cells
      } //for all contour values
    } //using scalar tree

  vtkDebugMacro(<<"Created: " 
               << newPts->GetNumberOfPoints() << " points, " 
               << newVerts->GetNumberOfCells() << " verts, " 
               << newLines->GetNumberOfCells() << " lines, " 
               << newPolys->GetNumberOfCells() << " triangles");
  //
  // Update ourselves.  Because we don't know up front how many verts, lines,
  // polys we've created, take care to reclaim memory. 
  //
  output->SetPoints(newPts);
  newPts->Delete();
  cellScalars->Delete();
  
  if (newVerts->GetNumberOfCells())
    {
    output->SetVerts(newVerts);
    }
  newVerts->Delete();

  if (newLines->GetNumberOfCells())
    {
    output->SetLines(newLines);
    }
  newLines->Delete();

  if (newPolys->GetNumberOfCells())
    {
    output->SetPolys(newPolys);
    }
  newPolys->Delete();

  this->Locator->Initialize();//releases leftover memory
  output->Squeeze();
}

//
// Special method handles structured points
//
#ifndef VTK_USE_PATENTED  
void vtkContourFilter::StructuredPointsContour(int vtkNotUsed(dim)) { }
#else
void vtkContourFilter::StructuredPointsContour(int dim)
{
  vtkPolyData *output;
  vtkPolyData *thisOutput = (vtkPolyData *)this->Output;
  int numContours=this->ContourValues->GetNumberOfContours();
  float *values=this->ContourValues->GetValues();

  if ( dim == 2 ) //marching squares
    {
    vtkMarchingSquares *msquares;
    int i;
    
    msquares = vtkMarchingSquares::New();
    msquares->SetInput((vtkStructuredPoints *)this->Input);
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
    mcubes->SetInput((vtkStructuredPoints *)this->Input);
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
#endif

// Specify a spatial locator for merging points. By default, 
// an instance of vtkMergePoints is used.
void vtkContourFilter::SetLocator(vtkPointLocator *locator)
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

void vtkContourFilter::CreateDefaultLocator()
{
  if ( this->Locator == NULL )
    {
    this->Locator = vtkMergePoints::New();
    }
}

void vtkContourFilter::PrintSelf(ostream& os, vtkIndent indent)
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

