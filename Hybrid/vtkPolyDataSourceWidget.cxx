/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkPolyDataSourceWidget.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkPolyDataSourceWidget.h"

#include "vtkDataSet.h"
#include "vtkProp3D.h"


vtkCxxRevisionMacro(vtkPolyDataSourceWidget, "1.4");

vtkPolyDataSourceWidget::vtkPolyDataSourceWidget() : vtk3DWidget()
{
  // child classes should call this constructor so that the vtk3DWidget()
  // constructor can set up some pertinent variables (e.g. Input and Prop3D)
}

void vtkPolyDataSourceWidget::PlaceWidget()
{
  float bounds[6];

  if ( this->Prop3D )
    {
    this->Prop3D->GetBounds(bounds);
    }
  else if ( this->Input )
    {
    this->Input->Update();
    // TODO: cleanup
    double *dbounds = this->Input->GetBounds();
    bounds[0] = (float)dbounds[0];
    bounds[1] = (float)dbounds[1];
    bounds[2] = (float)dbounds[2];
    bounds[3] = (float)dbounds[3];
    bounds[4] = (float)dbounds[4];
    bounds[5] = (float)dbounds[5];
    }
  else
    {
    // if Prop3D and Input aren't set, we assume that we're going to
    // look at what the user has already done with our polydata (and this
    // should happen in the child PlaceWidget(bounds), but we have to setup
    // some defaults for misbehaving child classes
    bounds[0] = -1.0;
    bounds[1] = 1.0;
    bounds[2] = -1.0;
    bounds[3] = 1.0;
    bounds[4] = -1.0;
    bounds[5] = 1.0;
    }
    
  this->PlaceWidget(bounds);
}

void vtkPolyDataSourceWidget::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}
