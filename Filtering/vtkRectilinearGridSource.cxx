/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkRectilinearGridSource.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkRectilinearGridSource.h"

#include "vtkObjectFactory.h"
#include "vtkRectilinearGrid.h"

vtkCxxRevisionMacro(vtkRectilinearGridSource, "1.18");

//----------------------------------------------------------------------------
vtkRectilinearGridSource::vtkRectilinearGridSource()
{
  this->vtkSource::SetNthOutput(0,vtkRectilinearGrid::New());
  // Releasing data for pipeline parallism.
  // Filters will know it is empty. 
  this->Outputs[0]->ReleaseData();
  this->Outputs[0]->Delete();
}

//----------------------------------------------------------------------------
vtkRectilinearGrid *vtkRectilinearGridSource::GetOutput(int idx)
{
  return (vtkRectilinearGrid *) this->vtkSource::GetOutput(idx); 
}

//----------------------------------------------------------------------------
vtkRectilinearGrid *vtkRectilinearGridSource::GetOutput()
{
  if (this->NumberOfOutputs < 1)
    {
    return NULL;
    }
  
  return (vtkRectilinearGrid *)(this->Outputs[0]);
}

//----------------------------------------------------------------------------
void vtkRectilinearGridSource::SetOutput(vtkRectilinearGrid *output)
{
  this->vtkSource::SetNthOutput(0, output);
}


//----------------------------------------------------------------------------
void vtkRectilinearGridSource::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}
