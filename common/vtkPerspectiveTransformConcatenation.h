/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkPerspectiveTransformConcatenation.h
  Language:  C++
  Date:      $Date$
  Version:   $Revision$
  Thanks:    Thanks to David G. Gobbi who developed this class.

Copyright (c) 1993-2000 Ken Martin, Will Schroeder, Bill Lorensen 
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
ARE DISCLAIMED. IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=========================================================================*/
// .NAME vtkPerspectiveTransformConcatenation - concatenation of perspective 
// transformations
// .SECTION Description
// This class is obsolete, all of its functionality can be found
// in vtkProjectionTransform.

#ifndef __vtkPerspectiveTransformConcatenation_h
#define __vtkPerspectiveTransformConcatenation_h

#include "vtkProjectionTransform.h"

class VTK_EXPORT vtkPerspectiveTransformConcatenation : 
  public vtkProjectionTransform
{
public:
  static vtkPerspectiveTransformConcatenation *New();

  vtkTypeMacro(vtkPerspectiveTransformConcatenation,vtkPerspectiveTransform);
  void PrintSelf(ostream& os, vtkIndent indent);
  
  // Description:
  // Make another transform of the same type.
  vtkGeneralTransform *MakeTransform();

protected:
  vtkPerspectiveTransformConcatenation();
  ~vtkPerspectiveTransformConcatenation();
  vtkPerspectiveTransformConcatenation(
	  const vtkPerspectiveTransformConcatenation&) {};
  void operator=(const vtkPerspectiveTransformConcatenation&) {};
};

#endif





