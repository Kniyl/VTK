/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkMesaProjectedPolyDataRayBounder.cxx
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
// Make sure this is first, so any includes of gl.h can be stoped if needed
#define VTK_IMPLEMENT_MESA_CXX

#include "GL/gl_mangle.h"
#include "GL/gl.h"

#include <math.h>
#include "vtkMesaProjectedPolyDataRayBounder.h"
#include "vtkRenderWindow.h"
#include "vtkMesaProperty.h"
#include "vtkMesaCamera.h"
#include "vtkMesaLight.h"
#include "vtkCuller.h"


// make sure this file is included before the #define takes place
// so we don't get two vtkMesaProjectedPolyDataRayBounder classes defined.
#include "vtkOpenGLProjectedPolyDataRayBounder.h"
#include "vtkMesaProjectedPolyDataRayBounder.h"

// Make sure vtkMesaProjectedPolyDataRayBounder is a copy of vtkOpenGLProjectedPolyDataRayBounder
// with vtkOpenGLProjectedPolyDataRayBounder replaced with vtkMesaProjectedPolyDataRayBounder
#define vtkOpenGLProjectedPolyDataRayBounder vtkMesaProjectedPolyDataRayBounder
#include "vtkOpenGLProjectedPolyDataRayBounder.cxx"
#undef vtkOpenGLProjectedPolyDataRayBounder

//-----------------------------------------------------------------------------
vtkMesaProjectedPolyDataRayBounder* vtkMesaProjectedPolyDataRayBounder::New()
{
  // First try to create the object from the vtkObjectFactory
  vtkObject* ret = vtkObjectFactory::CreateInstance("vtkMesaProjectedPolyDataRayBounder");
  if(ret)
    {
    return (vtkMesaProjectedPolyDataRayBounder*)ret;
    }
  // If the factory was unable to create the object, then create it here.
  return new vtkMesaProjectedPolyDataRayBounder;
}
