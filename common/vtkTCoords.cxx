/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkTCoords.cxx
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
#include "vtkTCoords.h"
#include "vtkObjectFactory.h"



//------------------------------------------------------------------------------
vtkTCoords* vtkTCoords::New()
{
  // First try to create the object from the vtkObjectFactory
  vtkObject* ret = vtkObjectFactory::CreateInstance("vtkTCoords");
  if(ret)
    {
    return (vtkTCoords*)ret;
    }
  // If the factory was unable to create the object, then create it here.
  return new vtkTCoords;
}




vtkTCoords *vtkTCoords::New(int dataType, int numComp)
{
  vtkTCoords *res = vtkTCoords::New();
  res->SetDataType(dataType);
  res->SetNumberOfComponents(numComp);
  return res;
}

// Construct object with an initial data array of type float and initial
// texture dimension of 2.
vtkTCoords::vtkTCoords()
{
  this->Data->SetNumberOfComponents(2);
}

// Set the data for this object. The tuple dimension must be consistent with
// the object.
void vtkTCoords::SetData(vtkDataArray *data)
{
  if ( data != this->Data && data != NULL )
    {
    if (data->GetNumberOfComponents() > 3 )
      {
      vtkErrorMacro(<<"Tuple dimension for texture coordinates must be <= 3");
      return;
      }
    this->Data->UnRegister(this);
    this->Data = data;
    this->Data->Register(this);
    this->Modified();
    }
}

// Given a list of pt ids, return an array of texture coordinates.
void vtkTCoords::GetTCoords(vtkIdList *ptIds, vtkTCoords *tc)
{
  int num = ptIds->GetNumberOfIds();

  tc->SetNumberOfTCoords(num);
  for (int i=0; i<num; i++)
    {
    tc->SetTCoord(i,this->GetTCoord(ptIds->GetId(i)));
    }
}

void vtkTCoords::PrintSelf(ostream& os, vtkIndent indent)
{
  vtkAttributeData::PrintSelf(os,indent);

  os << indent << "Number Of Texture Coordinates: " << this->GetNumberOfTCoords() << "\n";
  os << indent << "Number Of Texture Components: " << this->GetNumberOfComponents() << "\n";
}
