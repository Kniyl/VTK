/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkKitwareObjectFactory.h
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
// .NAME vtkKitwareObjectFactory - Object Factory for Kitware patented objects.
// .SECTION Description
// This is an object factory used to create Kitware patented objects.
// There is a KitwareFactory.dsp and KitwareFactory.dsw file to create
// the factory dll with the microsoft compiler.  Once the Factory is 
// created, put the resulting dll in VTK_AUTOLOAD_PATH.  
//
// .SECTION See Also
// vtkObjectFactory

#ifndef __vtkKitwareObjectFactory_h
#define __vtkKitwareObjectFactory_h

#include "vtkObjectFactory.h"

class VTK_PATENTED_EXPORT vtkKitwareObjectFactory : public vtkObjectFactory
{
public:
  static vtkKitwareObjectFactory *New() {return new vtkKitwareObjectFactory;};
  vtkTypeRevisionMacro(vtkKitwareObjectFactory,vtkObjectFactory);
  void PrintSelf(ostream& os, vtkIndent indent);  
  virtual const char* GetVTKSourceVersion();
protected:
  vtkKitwareObjectFactory() {};
  virtual vtkObject* CreateObject(const char* vtkclassname );
private:
  vtkKitwareObjectFactory(const vtkKitwareObjectFactory&);  // Not implemented.
  void operator=(const vtkKitwareObjectFactory&);  // Not implemented.
};

extern "C" VTK_PATENTED_EXPORT vtkObjectFactory* vtkLoad();
#endif
