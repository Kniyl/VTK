/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkCylinder.h
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

  Copyright (c) 1993-2002 Ken Martin, Will Schroeder, Bill Lorensen 
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkCylinder - implicit function for a cylinder
// .SECTION Description
// vtkCylinder computes the implicit function and function gradient for a
// cylinder. vtkCylinder is a concrete implementation of vtkImplicitFunction.
// Cylinder is centered at Center and axes of rotation is along the
// y-axis. (Use the superclass' vtkImplicitFunction transformation matrix if
// necessary to reposition.)

// .SECTION Caveats
// The cylinder is infinite in extent. To truncate the cylinder use the 
// vtkImplicitBoolean in combination with clipping planes.


#ifndef __vtkCylinder_h
#define __vtkCylinder_h

#include "vtkImplicitFunction.h"

class VTK_FILTERING_EXPORT vtkCylinder : public vtkImplicitFunction
{
public:
  vtkTypeRevisionMacro(vtkCylinder,vtkImplicitFunction);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description
  // Construct cylinder radius of 0.5.
  static vtkCylinder *New();

  // Description
  // Evaluate cylinder equation F(x,y,z) = (x-x0)^2 + (z-z0)^2 - R^2.
  float EvaluateFunction(float x[3]);
  float EvaluateFunction(float x, float y, float z)
    {return this->vtkImplicitFunction::EvaluateFunction(x, y, z); } ;

  // Description
  // Evaluate cylinder function gradient.
  void EvaluateGradient(float x[3], float g[3]);

  // Description:
  // Set/Get cylinder radius.
  vtkSetMacro(Radius,float);
  vtkGetMacro(Radius,float);

  // Description:
  // Set/Get cylinder center
  vtkSetVector3Macro(Center,float);
  vtkGetVectorMacro(Center,float,3);
protected:
  vtkCylinder();
  ~vtkCylinder() {};

  float Radius;
  float Center[3];

private:
  vtkCylinder(const vtkCylinder&);  // Not implemented.
  void operator=(const vtkCylinder&);  // Not implemented.
};

#endif


