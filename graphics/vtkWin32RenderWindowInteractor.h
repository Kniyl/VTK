/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkWin32RenderWindowInteractor.h
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
// .NAME vtkWin32RenderWindowInteractor - provide an event driven interface 
// to the renderer
// .SECTION Description
// vtkWin32RenderWindowInteractor is a convenience object that provides event 
// event bindings to common graphics functions. For example, camera and actor
// functions such as zoom-in/zoom-out, azimuth, roll, and pan. It is one of
// the window system specific subclasses of vtkRenderWindowInteractor. Please
// see vtkRenderWindowInteractor documentation for event bindings.
//
// Win32RenderWindowInteractor has an additional "animate" mode, which
// currently does not perform any animation, and has been left in the source
// code for compatibility.  Here both lower case and upper case will work
//
//<PRE>
// Additional keyboard bindings:
//    a - animate
//</PRE>
//
// .SECTION see also
// vtkRenderWindowInteractor vtkWin32OpenGLRenderWindow



#ifndef __vtkWin32RenderWindowInteractor_h
#define __vtkWin32RenderWindowInteractor_h

#include <stdlib.h>
#include "vtkRenderWindowInteractor.h"

class VTK_EXPORT vtkWin32RenderWindowInteractor : public vtkRenderWindowInteractor
{
public:

  // Description:
  // Construct object so that light follows camera motion.
  vtkWin32RenderWindowInteractor();
  ~vtkWin32RenderWindowInteractor();
  static vtkWin32RenderWindowInteractor *New() {
    return new vtkWin32RenderWindowInteractor;};
  const char *GetClassName() {return "vtkWin32RenderWindowInteractor";};
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Initialize the even handler
  virtual void Initialize();

  // Description:
  // Enable/Disable interactions.  By default interactors are enabled when
  // initialized.  Initialize() must be called prior to enabling/disabling
  // interaction. These methods are used when a window/widget is being
  // shared by multiple renderers and interactors.  This allows a "modal"
  // display where one interactor is active when its data is to be displayed
  // and all other interactors associated with the widget are disabled
  // when their data is not displayed.
  virtual void Enable();
  virtual void Disable();

  // Description:
  // This will start up the event loop and never return. If you
  // call this method it will loop processing events until the
  // application is exited.
  virtual void Start();
  void UpdateSize(int,int);

  // Description:
  // Provide implementaitons of the methods defined in
  // vtkRenderWindowInteractor. Generally the application developer should
  // not invoke these methods directly.
  virtual void StartRotate();
  virtual void EndRotate();
  virtual void StartZoom();
  virtual void EndZoom();
  virtual void StartPan();
  virtual void EndPan();
  virtual void StartSpin();
  virtual void EndSpin();
  virtual void StartDolly();
  virtual void EndDolly();
  virtual void StartUniformScale();
  virtual void EndUniformScale();

  virtual void StartAnimation();
  virtual void EndAnimation();

  //BTX
  friend VTK_EXPORT LRESULT CALLBACK vtkHandleMessage(HWND hwnd,UINT uMsg,
					   WPARAM w, LPARAM l);
  //ETX

  // Description:
  // Methods to set the default exit method for the class. This method is
  // only used if no instance level ExitMethod has been defined.  It is
  // provided as a means to control how an interactor is exited given
  // the various language bindings (tcl, Win32, etc.).
  static void SetClassExitMethod(void (*f)(void *), void *arg);
  static void SetClassExitMethodArgDelete(void (*f)(void *));
  
protected:
  HWND WindowId;
  UINT TimerId;
  WNDPROC OldProc;
  LPARAM LastPosition;
  
  //BTX
  // Description:
  // Class variables so an exit method can be defined for this class
  // (used to set different exit methods for various language bindings,
  // i.e. tcl, java, Win32)
  static void (*ClassExitMethod)(void *);
  static void (*ClassExitMethodArgDelete)(void *);
  static void *ClassExitMethodArg;
  //ETX
};

#endif


