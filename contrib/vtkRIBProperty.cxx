/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkRIBProperty.cxx
  Language:  C++
  Date:      $Date$
  Version:   $Revision$


Copyright (c) 1993-2001 Ken Martin, Will Schroeder, Bill Lorensen 
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.

 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

 * Neither name of Ken Martin, Will Schroeder, or Bill Lorensen nor the names
   of any contributors may be used to endorse or promote products derived
   from this software without specific prior written permission.

 * Modified source versions must be plainly marked as such, and must not be
   misrepresented as being the original software.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS ``AS IS''
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=========================================================================*/
#include "vtkRIBProperty.h"
#include "vtkObjectFactory.h"



//------------------------------------------------------------------------------
vtkRIBProperty* vtkRIBProperty::New()
{
  // First try to create the object from the vtkObjectFactory
  vtkObject* ret = vtkObjectFactory::CreateInstance("vtkRIBProperty");
  if(ret)
    {
    return (vtkRIBProperty*)ret;
    }
  // If the factory was unable to create the object, then create it here.
  return new vtkRIBProperty;
}




vtkRIBProperty::vtkRIBProperty ()
{
  this->Declarations = NULL;
  this->Parameters = NULL;
  this->SurfaceShader = new char[strlen("plastic") + 1];
  strcpy (this->SurfaceShader, "plastic");
  this->DisplacementShader = NULL;
  // create a vtkProperty that can be rendered
  this->Property = vtkProperty::New ();
}

vtkRIBProperty::~vtkRIBProperty()
{
  if (this->SurfaceShader)
    {
    delete [] this->SurfaceShader;
    }
  if (this->DisplacementShader)
    {
    delete [] this->DisplacementShader;
    }
  if (this->Declarations)
    {
    delete [] this->Declarations;
    }
  if (this->Property)
    {
    this->Property->Delete ();
    }
  if (this->Parameters)
    {
    delete [] this->Parameters;
    }
}

void vtkRIBProperty::Render(vtkActor *anActor, vtkRenderer *ren)
{
  int ref;
  
  // Copy this property's ivars into the property to be rendered
  ref = this->Property->GetReferenceCount();
  this->Property->DeepCopy(this);
  //this->Property->SetDeleteMethod(NULL);
  this->Property->SetReferenceCount(ref);
  
  // Render the property
  this->Property->Render (anActor, ren);
}

void vtkRIBProperty::SetVariable (char *variable, char *value)
{
  if (this->Declarations)
    {
    delete [] this->Declarations;
    }

  // format of line is: Declare "variable" "type"\n
  this->Declarations = new char [strlen ("Declare ") +
	                      strlen (variable) +
			      strlen (value) + 
			      8];

  sprintf (this->Declarations, "Declare \"%s\" \"%s\"\n", variable, value);
  this->Modified ();
}

void vtkRIBProperty::AddVariable (char *variable, char *value)
{
  if (this->Declarations == NULL)
    {
    this->SetVariable (variable, value);
    }
  else
    {
    char *newVariable = new char [strlen ("Declare ") +
	                          strlen (variable) +
		   	          strlen (value) + 
			          8];

    sprintf (newVariable, "Declare \"%s\" \"%s\"\n", variable, value);
    char *oldDeclarations = this->Declarations;

    this->Declarations = new char [strlen (oldDeclarations) + strlen (newVariable) + 1];
    strcpy (this->Declarations, oldDeclarations);
    strcat (this->Declarations, newVariable);
    delete [] oldDeclarations;
    delete [] newVariable;
    this->Modified ();
    }
}

void vtkRIBProperty::SetParameter (char *parameter, char *value)
{
  if (this->Parameters)
    {
    delete [] this->Parameters;
    }

  // format of line is: "parameter" "value"
  this->Parameters = new char [strlen (parameter) +
			      strlen (value) + 
			      7];

  sprintf (this->Parameters, " \"%s\" [%s]", parameter, value);
  this->Modified ();
}

void vtkRIBProperty::AddParameter (char *Parameter, char *value)
{
  if (this->Parameters == NULL)
    {
    this->SetParameter (Parameter, value);
    }
  else
    {
    char *newParameter = new char [strlen (Parameter) +
		   	          strlen (value) + 
			          7];

    sprintf (newParameter, " \"%s\" [%s]", Parameter, value);
    char *oldParameters = this->Parameters;

    this->Parameters = new char [strlen (oldParameters) + strlen (newParameter) + 1];
    strcpy (this->Parameters, oldParameters);
    strcat (this->Parameters, newParameter);
    delete [] oldParameters;
    delete [] newParameter;
    this->Modified ();
    }
}

char *vtkRIBProperty::GetParameters ()
{
  return this->Parameters;
}

char *vtkRIBProperty::GetDeclarations ()
{
  return this->Declarations;
}

void vtkRIBProperty::PrintSelf(ostream& os, vtkIndent indent)
{
  vtkProperty::PrintSelf(os,indent);
 
  if (this->SurfaceShader)
    {
    os << indent << "SurfaceShader: " << this->SurfaceShader << "\n";
    }
  else
    {
    os << indent << "SurfaceShader: (none)\n";
    }
  if (this->DisplacementShader)
    {
    os << indent << "DisplacementShader: " << this->DisplacementShader << "\n";
    }
  else
    {
    os << indent << "DisplacementShader: (none)\n";
    }
  if (this->Declarations)
    {
    os << indent << "Declarations: " << this->Declarations;
    }
  else
    {
    os << indent << "Declarations: (none)\n";
    }
  if (this->Parameters)
    {
    os << indent << "Parameters: " << this->Parameters;
    }
  else
    {
    os << indent << "Parameters: (none)\n";
    }

}

