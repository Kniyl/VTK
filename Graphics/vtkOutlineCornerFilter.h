/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkOutlineCornerFilter.h
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
// .NAME vtkOutlineCornerFilter - create wireframe outline corners for arbitrary data set
// .SECTION Description
// vtkOutlineCornerFilter is a filter that generates wireframe outline corners of any 
// data set. The outline consists of the eight corners of the dataset 
// bounding box.

#ifndef __vtkOutlineCornerFilter_h
#define __vtkOutlineCornerFilter_h

#include "vtkDataSetToPolyDataFilter.h"
class vtkOutlineCornerSource;

class VTK_GRAPHICS_EXPORT vtkOutlineCornerFilter : public vtkDataSetToPolyDataFilter
{
public:
  vtkTypeRevisionMacro(vtkOutlineCornerFilter,vtkDataSetToPolyDataFilter);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Construct outline corner filter with default corner factor = 0.2
  static vtkOutlineCornerFilter *New();

  // Description:
  // Set/Get the factor that controls the relative size of the corners
  // to the length of the corresponding bounds
  vtkSetClampMacro(CornerFactor, float, 0.001, 0.5);
  vtkGetMacro(CornerFactor, float);

protected:
  vtkOutlineCornerFilter();
  ~vtkOutlineCornerFilter();

  vtkOutlineCornerSource *OutlineCornerSource;
  void Execute();
  void ExecuteInformation();

  float CornerFactor;
private:
  vtkOutlineCornerFilter(const vtkOutlineCornerFilter&);  // Not implemented.
  void operator=(const vtkOutlineCornerFilter&);  // Not implemented.
};

#endif


