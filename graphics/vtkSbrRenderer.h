/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkSbrRenderer.h
  Language:  C++
  Date:      $Date$
  Version:   $Revision$


Copyright (c) 1993-1996 Ken Martin, Will Schroeder, Bill Lorensen.

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
// .NAME vtkSbrRenderer - HP starbase renderer
// .SECTION Description
// vtkSbrRenderer is a concrete implementation of the abstract class
// vtkRenderer. vtkSbrRenderer interfaces to the Hewlett-Packard starbase
// graphics library.

#ifndef __vtkSbrRenderer_h
#define __vtkSbrRenderer_h

#include <stdlib.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include "vtkRenderer.h"
#include "starbase.c.h"

class VTK_EXPORT vtkSbrRenderer : public vtkRenderer
{
public:
  vtkSbrRenderer();
  vtkSbrRenderer *New() {return new vtkSbrRenderer;};
  char *GetClassName() {return "vtkSbrRenderer";};
  void PrintSelf(ostream& os, vtkIndent indent);

  void Render(void);
  void ClearLights(void);

  int UpdateActors(void);
  int UpdateCameras(void);
  int UpdateLights(void);
  int UpdateVolumes(void);

  vtkGetMacro(Fd,int);
  vtkGetMacro(LightSwitch,int);
  vtkSetMacro(LightSwitch,int);
  virtual float *GetCenter();
  virtual void DisplayToView(); 
  virtual void ViewToDisplay(); 
  virtual int  IsInViewport(int x,int y); 

protected:
  int NumberOfLightsBound;
  int LightSwitch;
  int Fd;

};

#endif
