/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkMesaPolyDataMapper2D.h
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
// .NAME vtkWin32PolyDataMapper2D - 2D PolyData support for windows
// .SECTION Description
// vtkWin32PolyDataMapper2D provides 2D PolyData annotation support for 
// vtk under windows.  Normally the user should use vtkPolyDataMapper2D 
// which in turn will use this class.

// .SECTION See Also
// vtkPolyDataMapper2D

#ifndef __vtkMesaPolyDataMapper2D_h
#define __vtkMesaPolyDataMapper2D_h

#include "vtkPolyDataMapper2D.h"

class VTK_RENDERING_EXPORT vtkMesaPolyDataMapper2D : public vtkPolyDataMapper2D
{
public:
  vtkTypeRevisionMacro(vtkMesaPolyDataMapper2D,vtkPolyDataMapper2D);
  static vtkMesaPolyDataMapper2D *New();

  // Description:
  // Actually draw the poly data.
  void RenderOpaqueGeometry(vtkViewport* viewport, vtkActor2D* actor);

protected:
  vtkMesaPolyDataMapper2D() {};
  ~vtkMesaPolyDataMapper2D() {};
  
private:
  vtkMesaPolyDataMapper2D(const vtkMesaPolyDataMapper2D&);  // Not implemented.
  void operator=(const vtkMesaPolyDataMapper2D&);  // Not implemented.
};


#endif

