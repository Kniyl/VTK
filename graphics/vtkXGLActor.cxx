/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkXGLActor.cxx
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
#include "vtkXGLRenderer.h"
#include "vtkXGLActor.h"
#include "vtkMatrix4x4.h"

// Implement base class method.
void vtkXGLActor::Render(vtkRenderer *ren, vtkMapper *mapper)
{
  vtkMatrix4x4 *matrix = vtkMatrix4x4::New();
  Xgl_trans model_trans;

  // build transformation 
  this->GetMatrix(matrix);
  matrix->Transpose();

  // insert model transformation 
  xgl_object_get(*(((vtkXGLRenderer *)ren)->GetContext()),
                 XGL_CTX_GLOBAL_MODEL_TRANS, &model_trans);
  xgl_object_set(model_trans,  XGL_TRANS_DATA_TYPE, XGL_DATA_DBL, 0);
  xgl_transform_write(model_trans,(double (*)[4])(matrix->Element[0]));

  // send a render to the mapper; update pipeline
  mapper->Render(ren,this);
  matrix->Delete();
}

