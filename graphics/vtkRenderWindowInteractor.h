/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkRenderWindowInteractor.h
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
// .NAME vtkRenderWindowInteractor - base class for platform dependent
// implementations which handle routing of mouse/key/timer messages to
// vtkInterActorStyle and subclasses which handle all motion control
//
// .SECTION Description
// vtkRenderWindowInteractor has changed from previous implementations and
// now serves only as a shell to hold user preferences and route messages
// to vtkInterActorStyle. Callbacks are available for many Events.
// Platform specific subclasses should provide methods for
// CreateTimer/DestroyTimer, TerminateApp, and an event loop if required
// via Initialize/Start/Enable/Disable.

#ifndef __vtkRenderWindowInteractor_h
#define __vtkRenderWindowInteractor_h

#include "vtkObject.h"
#include "vtkRenderWindow.h"

// Timer flags for win32/X compatibility
#define VTKI_TIMER_FIRST  0
#define VTKI_TIMER_UPDATE 1

class vtkInteractorStyle;
class vtkPicker;
class vtkCellPicker;

class VTK_EXPORT vtkRenderWindowInteractor : public vtkObject
{
public:
  vtkRenderWindowInteractor();
  ~vtkRenderWindowInteractor();
  static vtkRenderWindowInteractor *New();
  const char *GetClassName() {return "vtkRenderWindowInteractor";};
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Prepare for handling events. This must be called before the
  // interactor will work.
  virtual void Initialize() {this->Initialized=1; this->Enable();
                             this->RenderWindow->Render();}

  // Description:
  // This Method detects loops of RenderWindow-Interactor,
  // so objects are freed properly.
  void UnRegister(vtkObject *o);

  // Description:
  // Start the event loop. This is provided so that you do not have to
  // implement your own event loop. You still can use your own
  // event loop if you want. Initialize should be called before Start.
  virtual void Start() {};

  // Description:
  // Enable/Disable interactions.  By default interactors are enabled when
  // initialized.  Initialize() must be called prior to enabling/disabling
  // interaction. These methods are used when a window/widget is being
  // shared by multiple renderers and interactors.  This allows a "modal"
  // display where one interactor is active when its data is to be displayed
  // and all other interactors associated with the widget are disabled
  // when their data is not displayed.
  virtual void Enable() { this->Enabled = 1; this->Modified();};
  virtual void Disable() { this->Enabled = 0; this->Modified();};
  vtkGetMacro(Enabled, int);

  // Description:
  // Set/Get the rendering window being controlled by this object.
  void SetRenderWindow(vtkRenderWindow *aren);
  vtkGetObjectMacro(RenderWindow,vtkRenderWindow);

  // Description:
  // Event loop notification member for Window size change
  virtual void UpdateSize(int x,int y);

  // Description:
  // This methods sets the Size ivar of the interactor without
  // actually changing the size of the window. Normally
  // application programmers would use UpdateSize if anything.
  // This is useful for letting someone else change the size of
  // the rendering window and just letting the interactor
  // know about the change.
  vtkSetVector2Macro(Size,int);
  vtkGetVector2Macro(Size,int);

  // Description:
  // Timer methods must be overridden by platform dependent subclasses.
  // flag is passed to indicate if this is first timer set or an update
  // as Win32 uses repeating timers, whereas X uses One shot more timer
  // if flag==VTKXI_TIMER_FIRST Win32 and X should createtimer
  // otherwise Win32 should exit and X should perform AddTimeOut()
  virtual bool CreateTimer(int timerflag)  { return true; };
  virtual bool DestroyTimer()              { return true; };

  // Description:
  // This function is called on 'q','e' keypress if exitmethod is not
  // specified and should be overidden by platform dependent subclasses
  // to provide a termination procedure if one is required.
  virtual void TerminateApp(void) {};

  // Description:
  // External switching between joystick/trackball/new? modes.
  virtual void SetInteractorStyle(vtkInteractorStyle *);
  vtkGetObjectMacro(InteractorStyle,vtkInteractorStyle);

  // Description:
  // Turn on/off the automatic repositioning of lights as the camera moves.
  vtkSetMacro(LightFollowCamera,int);
  vtkGetMacro(LightFollowCamera,int);
  vtkBooleanMacro(LightFollowCamera,int);

  // Description:
  // This method can be used by user callbacks to get the
  // x, y, coordinates of the current event.
  vtkSetVector2Macro(EventPosition,int);
  vtkGetVectorMacro(EventPosition,int,2);

  // Description:
  // Set/Get the desired update rate. This is used by vtkLODActor's to tell
  // them how quickly they need to render.  This update is in effect only
  // when the camera is being rotated, or zoomed.  When the interactor is
  // still, the StillUpdateRate is used instead.
  vtkSetClampMacro(DesiredUpdateRate,float,0.0001,VTK_LARGE_FLOAT);
  vtkGetMacro(DesiredUpdateRate,float);

  // Description:
  // Set/Get the desired update rate when movement has stopped.
  // See the SetDesiredUpdateRate method.
  vtkSetClampMacro(StillUpdateRate,float,0.0001,VTK_LARGE_FLOAT);
  vtkGetMacro(StillUpdateRate,float);

  // Description:
  // See whether interactor has been initialized yet.
  vtkGetMacro(Initialized,int);

  // Description:
  // Set the object used to perform pick operations. You can use this to
  // control what type of data is picked.
  void SetPicker(vtkPicker *picker);

  // Description:
  // Get the object used to perform pick operations.
  vtkGetObjectMacro(Picker,vtkPicker);

  // Description:
  // Create default picker. Used to create one when none is specified.
  virtual vtkPicker *CreateDefaultPicker();

  // Description:
  // Specify a method to be executed prior to the pick operation.
  void SetStartPickMethod(void (*f)(void *), void *arg);

  // Description:
  // Called when a void* argument is being discarded.  Lets the user free it.
  void SetStartPickMethodArgDelete(void (*f)(void *));

  // Description:
  // Specify a method to be executed after the pick operation.
  void SetEndPickMethod(void (*f)(void *), void *arg);

  // Description:
  // Called when a void* argument is being discarded.  Lets the user free it.
  void SetEndPickMethodArgDelete(void (*f)(void *));

  // Description:
  // Set the user method. This method is invoked on a "u" keypress.
  void SetUserMethod(void (*f)(void *), void *arg);

  // Description:
  // Called when a void* argument is being discarded.  Lets the user free it.
  void SetUserMethodArgDelete(void (*f)(void *));

  // Description:
  // Set the exit method. This method is invoked on a "e" or "q" keypress.
  void SetExitMethod(void (*f)(void *), void *arg);

  // Description:
  // Called when a void* argument is being discarded.  Lets the user free it.
  void SetExitMethodArgDelete(void (*f)(void *));

  // Description:
  // These methods correspond to the the Exit, User and Pick
  // callbacks. They allow for the Style to invoke them.
  virtual void ExitCallback();
  virtual void UserCallback();
  virtual void StartPickCallback();
  virtual void EndPickCallback();
  
  // Description:
  // Render the scene. Just pass the render call on to the RenderWindow.
  void Render();
  
protected:
  vtkRenderWindow    *RenderWindow;
  vtkInteractorStyle *InteractorStyle;
  // used to track picked objects in actor mode
  // reason for existence: user may use any kind of picker.  Interactor
  //    need the high precision of cell picker at all time.
  vtkPicker          *Picker;

  //
  int   Initialized;
  int   Enabled;
  int   Style;
  int   LightFollowCamera;
  int   ActorMode;
  float DesiredUpdateRate;
  float StillUpdateRate;
  int   EventPosition[2], Size[2];
  
  // user methods that can be used to override default behaviour
  void (*StartPickMethod)(void *);
  void (*StartPickMethodArgDelete)(void *);
  void *StartPickMethodArg;
  void (*EndPickMethod)(void *);
  void (*EndPickMethodArgDelete)(void *);
  void *EndPickMethodArg;

  void (*UserMethod)(void *);
  void (*UserMethodArgDelete)(void *);
  void *UserMethodArg;

  void (*ExitMethod)(void *);
  void (*ExitMethodArgDelete)(void *);
  void *ExitMethodArg;
};

#endif
