/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkLinearTransformInverse.cxx
  Language:  C++
  Date:      $Date$
  Version:   $Revision$
  Thanks:    Thanks to David G. Gobbi who developed this class.

Copyright (c) 1993-2000 Ken Martin, Will Schroeder, Bill Lorensen.

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

#include "vtkLinearTransformInverse.h"
#include "vtkMath.h"
#include "vtkObjectFactory.h"

//----------------------------------------------------------------------------
vtkLinearTransformInverse* vtkLinearTransformInverse::New()
{
  // First try to create the object from the vtkObjectFactory
  vtkObject* ret = vtkObjectFactory::CreateInstance("vtkLinearTransformInverse");
  if(ret)
    {
    return (vtkLinearTransformInverse*)ret;
    }
  // If the factory was unable to create the object, then create it here.
  return new vtkLinearTransformInverse;
}

//----------------------------------------------------------------------------
vtkLinearTransformInverse::vtkLinearTransformInverse()
{
  this->TransformType = VTK_INVERSE_TRANSFORM | VTK_LINEAR_TRANSFORM;

  this->Transform = NULL;
  this->UpdateRequired = 0;
}

//----------------------------------------------------------------------------
vtkLinearTransformInverse::~vtkLinearTransformInverse()
{
  if (this->Transform)
    {
    this->Transform->Delete();
    }
}

//----------------------------------------------------------------------------
void vtkLinearTransformInverse::PrintSelf(ostream& os, vtkIndent indent)
{
  vtkLinearTransform::PrintSelf(os,indent);

  os << indent << "Transform: " << this->Transform << "\n";
  if (this->Transform)
    {
    this->Transform->PrintSelf(os,indent.GetNextIndent());
    }
}

//----------------------------------------------------------------------------
void vtkLinearTransformInverse::SetInverse(vtkLinearTransform *trans)
{
  if (this == trans)
    {
    vtkErrorMacro(<<"SetInverse: A transform cannot be its own inverse!");
    return;
    }
  if (this->MyInverse == trans)
    {
    return;
    }
  if (this->MyInverse)
    {
    this->MyInverse->Delete();
    this->Transform->Delete();
    }
  this->MyInverse = trans;
  trans->Register(this);
  this->Transform = (vtkLinearTransform *)trans->MakeTransform();
  this->UpdateRequired = 1;
  this->Modified();
}

//----------------------------------------------------------------------------
vtkGeneralTransform *vtkLinearTransformInverse::GetInverse()
{
  return (vtkLinearTransform *)this->MyInverse;
}

//----------------------------------------------------------------------------
vtkLinearTransform *vtkLinearTransformInverse::GetTransform()
{
  return this->Transform;
}

//----------------------------------------------------------------------------
void vtkLinearTransformInverse::Identity()
{
  if (this->MyInverse == NULL)
    {
    vtkErrorMacro(<< "Identity: Inverse has not been set");
    return;
    }
  this->MyInverse->Identity();
}

//----------------------------------------------------------------------------
void vtkLinearTransformInverse::Inverse()
{
  if (this->MyInverse == NULL)
    {
    vtkErrorMacro(<< "Inverse: Inverse has not been set");
    return;
    }
  this->MyInverse->Inverse();
}

//----------------------------------------------------------------------------
vtkGeneralTransform *vtkLinearTransformInverse::MakeTransform()
{
  if (this->MyInverse == NULL)
    {
    vtkErrorMacro(<< "MakeTransform: Inverse has not been set");
    return NULL;
    }
  return this->MyInverse->MakeTransform();
}

//----------------------------------------------------------------------------
void vtkLinearTransformInverse::DeepCopy(vtkGeneralTransform *transform)
{
  if (this->MyInverse == NULL)
    {
    vtkErrorMacro(<< "DeepCopy: Inverse has not been set");
    return;
    }
  this->MyInverse->DeepCopy(transform);
  this->MyInverse->Inverse();
}

//----------------------------------------------------------------------------
void vtkLinearTransformInverse::Update()
{
  if (this->MyInverse == NULL)
    {
    vtkErrorMacro(<< "Update: Inverse has not been set");
    return;
    }

  this->MyInverse->Update();

  if (this->MyInverse->GetMTime() > 
      this->Transform->GetMTime() || this->UpdateRequired)
    {
    this->Transform->DeepCopy(this->MyInverse);
    this->Transform->Inverse();
    this->UpdateRequired = 0;
    }

  this->Transform->Update();
  this->Matrix->DeepCopy(this->Transform->GetMatrixPointer());
}

//----------------------------------------------------------------------------
unsigned long vtkLinearTransformInverse::GetMTime()
{
  unsigned long result = this->vtkLinearTransform::GetMTime();
  unsigned long mtime;

  if (this->MyInverse)
    {
    mtime = this->MyInverse->GetMTime();
    if (mtime > result)
      {
      result = mtime;
      }
    }
  return result;
}



