/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkImageInformation.cxx
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
#include "vtkImageInformation.h"
#include "vtkObjectFactory.h"



//------------------------------------------------------------------------------
vtkImageInformation* vtkImageInformation::New()
{
  // First try to create the object from the vtkObjectFactory
  vtkObject* ret = vtkObjectFactory::CreateInstance("vtkImageInformation");
  if(ret)
    {
    return (vtkImageInformation*)ret;
    }
  // If the factory was unable to create the object, then create it here.
  return new vtkImageInformation;
}




//----------------------------------------------------------------------------
// Construct a new vtkImageInformation 
vtkImageInformation::vtkImageInformation()
{
  this->Spacing[0] = 1.0;
  this->Spacing[1] = 1.0;
  this->Spacing[2] = 1.0;

  this->Origin[0] = 0.0;
  this->Origin[1] = 0.0;
  this->Origin[2] = 0.0;
  
  this->ScalarType = VTK_VOID;
  this->NumberOfScalarComponents = 1;
}


//----------------------------------------------------------------------------
void vtkImageInformation::PrintSelf(ostream& os, vtkIndent indent)
{
  vtkStructuredInformation::PrintSelf(os, indent);
  
  os << indent << "NumberOfScalarComponents: " << 
    this->NumberOfScalarComponents << endl;
  os << indent << "ScalarType: " << this->ScalarType << endl;

  os << indent << "Spacing: (" << this->Spacing[0] << ", "
                               << this->Spacing[1] << ", "
                               << this->Spacing[2] << ")\n";

  os << indent << "Origin: (" << this->Origin[0] << ", "
                              << this->Origin[1] << ", "
                              << this->Origin[2] << ")\n";
}


//----------------------------------------------------------------------------
int vtkImageInformation::GetClassCheck(char *className)
{
  if (strcmp(className, "vtkImageInformation") == 0)
    {
    return 1;
    }
  // check superclass
  if (this->vtkStructuredInformation::GetClassCheck(className))
    {
    return 1;
    }
  
  return 0;
}

//----------------------------------------------------------------------------
void vtkImageInformation::Copy(vtkDataInformation *in)
{
  this->vtkStructuredInformation::Copy(in);
  
  if (in->GetClassCheck("vtkImageInformation"))
    {
    vtkImageInformation *info = (vtkImageInformation*)(in);
    this->SetOrigin(info->GetOrigin());
    this->SetSpacing(info->GetSpacing());
    this->SetScalarType(info->GetScalarType());
    this->SetNumberOfScalarComponents(info->GetNumberOfScalarComponents());
    }
}

//----------------------------------------------------------------------------
void vtkImageInformation::WriteSelf(ostream& os)
{
  this->vtkStructuredInformation::WriteSelf(os);

  os << this->Spacing[0] << " " << this->Spacing[1] << " " 
     << this->Spacing[2] << " ";
  os << this->Origin[0] << " " << this->Origin[1] << " " 
     << this->Origin[2] << " ";
  os << this->ScalarType << " " ;
  os << this->NumberOfScalarComponents << " " ;
}

//----------------------------------------------------------------------------
void vtkImageInformation::ReadSelf(istream& is)
{
  this->vtkStructuredInformation::ReadSelf(is);

  is >> this->Spacing[0] >> this->Spacing[1] >> this->Spacing[2] ;
  is >> this->Origin[0] >> this->Origin[1] >> this->Origin[2] ;
  is >> this->ScalarType ;
  is >> this->NumberOfScalarComponents ;
}



  



