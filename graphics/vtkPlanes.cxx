/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkPlanes.cxx
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
#include "vtkPlanes.h"
#include "vtkPlane.h"
#include "vtkCamera.h"
#include "vtkObjectFactory.h"



//------------------------------------------------------------------------------
vtkPlanes* vtkPlanes::New()
{
  // First try to create the object from the vtkObjectFactory
  vtkObject* ret = vtkObjectFactory::CreateInstance("vtkPlanes");
  if(ret)
    {
    return (vtkPlanes*)ret;
    }
  // If the factory was unable to create the object, then create it here.
  return new vtkPlanes;
}




vtkPlanes::vtkPlanes()
{
  this->Points = NULL;
  this->Normals = NULL;
  this->Plane = vtkPlane::New();

  for (int i=0; i<24; i++)
    {
    this->Planes[i] = 0.0;
    }
}

vtkPlanes::~vtkPlanes()
{
  if ( this->Points )
    {
    this->Points->UnRegister(this);
    }
  if ( this->Normals )
    {
    this->Normals->UnRegister(this);
    }
  this->Plane->Delete();
}

// Evaluate plane equations. Return smallest absolute value.
float vtkPlanes::EvaluateFunction(float x[3])
{
  int numPlanes, i;
  float val, maxVal;

  if ( !this->Points || ! this->Normals )
    {
    vtkErrorMacro(<<"Please define points and/or normals!");
    return VTK_LARGE_FLOAT;
    }

  if ( (numPlanes=this->Points->GetNumberOfPoints()) != this->Normals->GetNumberOfNormals() )
    {
    vtkErrorMacro(<<"Number of normals/points inconsistent!");
    return VTK_LARGE_FLOAT;
    }

  for (maxVal=-VTK_LARGE_FLOAT, i=0; i < numPlanes; i++)
    {
    val = this->Plane->Evaluate(this->Normals->GetNormal(i),
			 this->Points->GetPoint(i), x);
    if (val > maxVal )
      {
      maxVal = val;
      }
    }

  return maxVal;
}

// Evaluate planes gradient.
void vtkPlanes::EvaluateGradient(float x[3], float n[3])
{
  int numPlanes, i;
  float val, maxVal, *nTemp;

  if ( !this->Points || ! this->Normals )
    {
    vtkErrorMacro(<<"Please define points and/or normals!");
    return;
    }

  if ( (numPlanes=this->Points->GetNumberOfPoints()) != this->Normals->GetNumberOfNormals() )
    {
    vtkErrorMacro(<<"Number of normals/points inconsistent!");
    return;
    }

  for (maxVal=-VTK_LARGE_FLOAT, i=0; i < numPlanes; i++)
    {
    nTemp = this->Normals->GetNormal(i);
    val = this->Plane->Evaluate(nTemp,this->Points->GetPoint(i), x);
    if ( val > maxVal )
      {
      maxVal = val;
      n[0] = nTemp[0];
      n[1] = nTemp[1];
      n[2] = nTemp[2];
      }
    }
}

void vtkPlanes::SetFrustumPlanes(float aspect, vtkCamera *camera)
{
  int i;
  float planes[24], *plane, n[3], x[3];
  
  // Get the planes and load them into the implicit function
  camera->GetFrustumPlanes(aspect,planes);
  for (i=0; i<24; i++)
    {
    if ( this->Planes[i] != planes[i] )
      {
      break;
      }
    }
  if ( i >= 24 )
    {
    return; //same as before don't modify
    }

  // okay, need to allocate stuff
  this->Modified();
  vtkPoints *pts = vtkPoints::New();
  vtkNormals *normals = vtkNormals::New();

  pts->SetNumberOfPoints(6);
  normals->SetNumberOfNormals(6);
  this->SetPoints(pts);
  this->SetNormals(normals);

  for (i=0; i<6; i++)
    {
    plane = planes + 4*i;
    n[0] = -plane[0];
    n[1] = -plane[1];
    n[2] = -plane[2];
    x[0] = x[1] = x[2] = 0.0;
    if ( n[0] != 0.0 )
      {
      x[0] = plane[3] / n[0];
      }
    else if ( n[1] != 0.0 )
      {
      x[1] = plane[3] / n[1];
      }
    else
      {
      x[2] = plane[3] / n[2];
      }
    pts->SetPoint(i,x);
    normals->SetNormal(i,n);
    }
  
  pts->Delete(); //ok reference counting
  normals->Delete();
}

void vtkPlanes::PrintSelf(ostream& os, vtkIndent indent)
{
  int numPlanes;

  vtkImplicitFunction::PrintSelf(os,indent);

  if ( this->Points && (numPlanes=this->Points->GetNumberOfPoints()) > 0 )
    {
    os << indent << "Number of Planes: " << numPlanes << "\n";
    }
  else
    {
    os << indent << "No Planes Defined.\n";
    }

  if ( this->Normals )
    {
    os << indent << "Normals: " << this->Normals << "\n";
    }
  else
    {
    os << indent << "Normals: (none)\n";
    }
}
