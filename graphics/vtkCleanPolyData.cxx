/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkCleanPolyData.cxx
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
#include "vtkCleanPolyData.h"
#include "vtkMergePoints.h"

// Construct object with initial tolerance of 0.0.
vtkCleanPolyData::vtkCleanPolyData()
{
  this->Tolerance = 0.0;
  this->Locator = NULL;
}

vtkCleanPolyData::~vtkCleanPolyData()
{
  if ( this->Locator )
    {
    this->Locator->UnRegister(this);
    this->Locator = NULL;
    }
}

void vtkCleanPolyData::Execute()
{
  vtkPolyData *input=(vtkPolyData *)this->Input;
  vtkPointData *pd=input->GetPointData();
  int numPts=input->GetNumberOfPoints();
  vtkPoints *inPts;
  vtkPoints *newPts;
  int numNewPts;
  int *updatedPts= new int[input->GetMaxCellSize()];
  int i, ptId;
  int npts, *pts;
  float *x;
  vtkCellArray *inVerts=input->GetVerts(), *newVerts=NULL;
  vtkCellArray *inLines=input->GetLines(), *newLines=NULL;
  vtkCellArray *inPolys=input->GetPolys(), *newPolys=NULL;
  vtkCellArray *inStrips=input->GetStrips(), *newStrips=NULL;
  vtkPolyData *output = this->GetOutput();
  vtkPointData *outputPD = output->GetPointData();
  
  vtkDebugMacro(<<"Cleaning data");

  if ( numPts < 1 || (inPts=input->GetPoints()) == NULL )
    {
    vtkErrorMacro(<<"No data to clean!");
    return;
    }

  outputPD->CopyAllocate(pd);
  this->CreateDefaultLocator();

  // Initialize; compute absolute tolerance from relative given
  this->Locator->SetTolerance(this->Tolerance*input->GetLength());
  newPts = vtkPoints::New();
  newPts->Allocate(numPts);
  this->Locator->InitPointInsertion (newPts, input->GetBounds());

  //
  // Begin to adjust topology.
  //
  // Vertices are renumbered and we remove duplicate vertices
  if ( inVerts->GetNumberOfCells() > 0 )
    {
    newVerts = vtkCellArray::New();
    newVerts->Allocate(inVerts->GetSize());

    for (inVerts->InitTraversal(); inVerts->GetNextCell(npts,pts); )
      {
      for ( numNewPts=0, i=0; i < npts; i++ )
        {
        x = inPts->GetPoint(pts[i]);

        if ( (ptId=this->Locator->IsInsertedPoint(x)) < 0 )
          {
          ptId = this->Locator->InsertNextPoint(x);
          updatedPts[numNewPts++] = ptId;
          outputPD->CopyData(pd,pts[i],ptId);
          }
        }
      if ( numNewPts > 0 )
	{
	newVerts->InsertNextCell(numNewPts,updatedPts);
	}
      }
    newVerts->Squeeze();
    }
  this->UpdateProgress(0.25);

  // lines reduced to one point are eliminated
  if ( inLines->GetNumberOfCells() > 0 )
    {
    newLines = vtkCellArray::New();
    newLines->Allocate(inLines->GetSize());

    for (inLines->InitTraversal(); inLines->GetNextCell(npts,pts); )
      {
      for ( numNewPts=0, i=0; i < npts; i++ )
        {
        x = inPts->GetPoint(pts[i]);

        if ( (ptId=this->Locator->IsInsertedPoint(x)) < 0 )
          {
          ptId = this->Locator->InsertNextPoint(x);
          outputPD->CopyData(pd,pts[i],ptId);
          }

        if ( i == 0 || ptId != updatedPts[numNewPts-1] )
          {
          updatedPts[numNewPts++] = ptId;
          }
        }

      if ( numNewPts > 1 )
	{
	newLines->InsertNextCell(numNewPts,updatedPts);
	}
      }
    newLines->Squeeze();
    vtkDebugMacro(<<"Removed " << inLines->GetNumberOfCells() -
                 newLines->GetNumberOfCells() << " lines");
    }
  this->UpdateProgress(0.50);

  // polygons reduced to two points or less are eliminated
  if ( inPolys->GetNumberOfCells() > 0 )
    {
    newPolys = vtkCellArray::New();
    newPolys->Allocate(inPolys->GetSize());
    for (inPolys->InitTraversal(); inPolys->GetNextCell(npts,pts); )
      {
      for ( numNewPts=0, i=0; i < npts; i++ )
        {
        x = inPts->GetPoint(pts[i]);

        if ( (ptId=this->Locator->IsInsertedPoint(x)) < 0 )
          {
          ptId = this->Locator->InsertNextPoint(x);
          outputPD->CopyData(pd,pts[i],ptId);
          }

	// check for duplicate points
        if ( i == 0 || ptId != updatedPts[numNewPts-1] )
          {
          updatedPts[numNewPts++] = ptId;
          }
        }//for points in polygon
      
      if ( numNewPts > 2 && updatedPts[0] == updatedPts[numNewPts-1] )
	{
	numNewPts--;
	}
      if ( numNewPts > 2 )
	{
	newPolys->InsertNextCell(numNewPts,updatedPts);
	}
      }

    newPolys->Squeeze();
    vtkDebugMacro(<<"Removed " << inPolys->GetNumberOfCells() -
                 newPolys->GetNumberOfCells() << " polys");
    }
  this->UpdateProgress(0.75);

  // triangle strips reduced to two points or less are eliminated
  if ( inStrips->GetNumberOfCells() > 0 ) 
    {
    newStrips = vtkCellArray::New();
    newStrips->Allocate(inStrips->GetSize());

    for (inStrips->InitTraversal(); inStrips->GetNextCell(npts,pts); )
      {
      for ( numNewPts=0, i=0; i < npts; i++ )
        {
        x = inPts->GetPoint(pts[i]);

        if ( (ptId=this->Locator->IsInsertedPoint(x)) < 0 )
          {
          ptId = this->Locator->InsertNextPoint(x);
          outputPD->CopyData(pd,pts[i],ptId);
          }

        if ( i == 0 || ptId != updatedPts[numNewPts-1] )
          {
          updatedPts[numNewPts++] = ptId;
          }
        }

      if ( numNewPts > 2 )
	{
	newStrips->InsertNextCell(numNewPts,updatedPts);
	}
      }

    newStrips->Squeeze();
    vtkDebugMacro(<<"Removed " << inStrips->GetNumberOfCells() -
                 newStrips->GetNumberOfCells() << " strips");
    }

  vtkDebugMacro(<<"Removed " << numPts - newPts->GetNumberOfPoints()
                << " points");

  // Update ourselves and release memory
  //
  delete [] updatedPts;
  this->Locator->Initialize(); //release memory.

  output->SetPoints(newPts);
  newPts->Delete();
  if (newVerts)
    {
    output->SetVerts(newVerts);
    newVerts->Delete();
    }
  if (newLines)
    {
    output->SetLines(newLines);
    newLines->Delete();
    }
  if (newPolys)
    {
    output->SetPolys(newPolys);
    newPolys->Delete();
    }
  if (newStrips)
    {
    output->SetStrips(newStrips);
    newStrips->Delete();
    }
}

// Specify a spatial locator for speeding the search process. By
// default an instance of vtkLocator is used.
void vtkCleanPolyData::SetLocator(vtkPointLocator *locator)
{
  if ( this->Locator == locator)
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

// Method manages creation of locators. It takes into account the potential
// change of tolerance (zero to non-zero).
void vtkCleanPolyData::CreateDefaultLocator()
{
  if ( this->Locator == NULL)
    {
    if ( this->Tolerance <= 0.0 )
      {
      this->Locator = vtkMergePoints::New();
      }
    else
      {
      this->Locator = vtkPointLocator::New();
      }
    }
  else
    {
    if ( !strcmp(this->Locator->GetClassName(),"vtkMergePoints") &&
	 this->Tolerance > 0.0 )
      {
      vtkWarningMacro(<<"Trying to merge points with non-zero tolerance using vtkMergePoints");
      }
    }
}

void vtkCleanPolyData::PrintSelf(ostream& os, vtkIndent indent)
{
  vtkPolyDataToPolyDataFilter::PrintSelf(os,indent);

  os << indent << "Tolerance: " << this->Tolerance << "\n";

  if ( this->Locator )
    {
    os << indent << "Locator: " << this->Locator << "\n";
    }
  else
    {
    os << indent << "Locator: (none)\n";
    }
}

unsigned long int vtkCleanPolyData::GetMTime()
{
  unsigned long mTime=this->vtkObject::GetMTime();
  unsigned long time;

  if ( this->Locator != NULL )
    {
    time = this->Locator->GetMTime();
    mTime = ( time > mTime ? time : mTime );
    }
  return mTime;
}
