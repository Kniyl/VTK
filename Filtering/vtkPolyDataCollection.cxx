/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkPolyDataCollection.cxx
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
#include "vtkPolyDataCollection.h"
#include "vtkObjectFactory.h"

vtkCxxRevisionMacro(vtkPolyDataCollection, "1.7");
vtkStandardNewMacro(vtkPolyDataCollection);

//----------------------------------------------------------------------------
void vtkPolyDataCollection::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}
