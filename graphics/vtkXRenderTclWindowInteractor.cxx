/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkXRenderTclWindowInteractor.cxx
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <X11/X.h>
#include <X11/keysym.h>
#include "vtkXRenderWindowInteractor.h"
#include "vtkXRenderWindow.h"
#include "vtkActor.h"
#include <X11/Shell.h>
#include <math.h>
#include "tk.h"

// steal the first two elements of the TkMainInfo stuct
// we don't care about the rest of the elements.
typedef struct TkMainInfo {
  int refCount;
  struct TkWindow *winPtr;
};

extern TkMainInfo *tkMainWindowList;

// returns 1 if done
static int vtkTclEventProc(XtPointer clientData,XEvent *event)
{
  Boolean ctd;
  vtkXRenderWindow *rw;
      
  rw = (vtkXRenderWindow *)
    (((vtkXRenderWindowInteractor *)clientData)->GetRenderWindow());
  
  if (rw->GetWindowId() == ((XAnyEvent *)event)->window)
    {
    vtkXRenderWindowInteractorCallback((Widget)NULL,clientData, event, &ctd);
    ctd = 0;
    }
  else
    {
    ctd = 1;
    }

  return !ctd;
}

static void vtkXTclTimerProc(ClientData clientData)
{
  XtIntervalId id;
  vtkXRenderWindowInteractorTimer((XtPointer)clientData,&id);
}

// states
#define VTKXI_START  0
#define VTKXI_ROTATE 1
#define VTKXI_ZOOM   2
#define VTKXI_PAN    3

// Description:
// Construct object so that light follows camera motion.
vtkXRenderWindowInteractor::vtkXRenderWindowInteractor()
{
  this->State = VTKXI_START;
  this->App = 0;
  this->top = 0;
}

vtkXRenderWindowInteractor::~vtkXRenderWindowInteractor()
{
  if (this->Initialized)
    {
    Tk_DeleteGenericHandler((Tk_GenericProc *)vtkTclEventProc,
			    (ClientData)this);
    }
}

void  vtkXRenderWindowInteractor::SetWidget(Widget foo)
{
  this->top = foo;
} 

void  vtkXRenderWindowInteractor::Start()
{
  Tk_MainLoop();
}

// Description:
// Initializes the event handlers
void vtkXRenderWindowInteractor::Initialize(XtAppContext app)
{
  this->App = app;

  this->Initialize();
}

// Description:
// Begin processing keyboard strokes.
void vtkXRenderWindowInteractor::Initialize()
{
  vtkXRenderWindow *ren;
  int depth;
  Colormap cmap;
  Visual  *vis;
  int *size;
  int *position;

  // make sure we have a RenderWindow and camera
  if ( ! this->RenderWindow)
    {
    vtkErrorMacro(<<"No renderer defined!");
    return;
    }

  this->Initialized = 1;
  ren = (vtkXRenderWindow *)(this->RenderWindow);

  // use the same display as tcl/tk
  ren->SetDisplayId(Tk_Display(tkMainWindowList->winPtr));
  this->DisplayId = ren->GetDisplayId();
  
  // get the info we need from the RenderingWindow
  depth   = ren->GetDesiredDepth();
  cmap    = ren->GetDesiredColormap();
  vis     = ren->GetDesiredVisual();
  size    = ren->GetSize();
  position= ren->GetPosition();
  
  size = ren->GetSize();
  ren->Render();
  this->WindowId = ren->GetWindowId();
  size = ren->GetSize();

  this->Size[0] = size[0];
  this->Size[1] = size[1];

  XSelectInput(this->DisplayId,this->WindowId,
	       KeyPressMask | ButtonPressMask | ExposureMask |
		    StructureNotifyMask | ButtonReleaseMask);

  // add in tcl init stuff
  Tk_CreateGenericHandler((Tk_GenericProc *)vtkTclEventProc,(ClientData)this);
}

void vtkXRenderWindowInteractor::PrintSelf(ostream& os, vtkIndent indent)
{
  vtkRenderWindowInteractor::PrintSelf(os,indent);
}


void  vtkXRenderWindowInteractor::UpdateSize(int x,int y)
{
  // if the size changed send this on to the RenderWindow
  if ((x != this->Size[0])||(y != this->Size[1]))
    {
    this->Size[0] = x;
    this->Size[1] = y;
    this->RenderWindow->SetSize(x,y);
    }

}
 
void  vtkXRenderWindowInteractor::StartRotate()
{
  if (this->State != VTKXI_START) return;
  this->State = VTKXI_ROTATE;

  this->RenderWindow->SetDesiredUpdateRate(this->DesiredUpdateRate);
  Tk_CreateTimerHandler(10,vtkXTclTimerProc,(ClientData)this);
}
void  vtkXRenderWindowInteractor::EndRotate()
{
  if (this->State != VTKXI_ROTATE) return;
  this->State = VTKXI_START;
  this->RenderWindow->SetDesiredUpdateRate(this->StillUpdateRate);
  this->RenderWindow->Render();
}

void  vtkXRenderWindowInteractor::StartZoom()
{
  if (this->State != VTKXI_START) return;
  this->State = VTKXI_ZOOM;
  this->RenderWindow->SetDesiredUpdateRate(this->DesiredUpdateRate);
  Tk_CreateTimerHandler(10,vtkXTclTimerProc,(ClientData)this);
}
void  vtkXRenderWindowInteractor::EndZoom()
{
  if (this->State != VTKXI_ZOOM) return;
  this->State = VTKXI_START;
  this->RenderWindow->SetDesiredUpdateRate(this->StillUpdateRate);
  this->RenderWindow->Render();
}

void  vtkXRenderWindowInteractor::StartPan()
{
  float *FocalPoint;
  float *Result;

  if (this->State != VTKXI_START) return;

  this->State = VTKXI_PAN;

  // calculate the focal depth since we'll be using it a lot
  FocalPoint = this->CurrentCamera->GetFocalPoint();
      
  this->CurrentRenderer->SetWorldPoint(FocalPoint[0],FocalPoint[1],
				       FocalPoint[2],1.0);
  this->CurrentRenderer->WorldToDisplay();
  Result = this->CurrentRenderer->GetDisplayPoint();
  this->FocalDepth = Result[2];

  this->RenderWindow->SetDesiredUpdateRate(this->DesiredUpdateRate);
  Tk_CreateTimerHandler(10,vtkXTclTimerProc,(ClientData)this);
}
void  vtkXRenderWindowInteractor::EndPan()
{
  if (this->State != VTKXI_PAN) return;
  this->State = VTKXI_START;
  this->RenderWindow->SetDesiredUpdateRate(this->StillUpdateRate);
  this->RenderWindow->Render();
}

void vtkXRenderWindowInteractorCallback(Widget vtkNotUsed(w),
					XtPointer client_data, 
					XEvent *event, 
					Boolean *vtkNotUsed(ctd))
{
  vtkXRenderWindowInteractor *me;
  XEvent marker;
  
  me = (vtkXRenderWindowInteractor *)client_data;

  switch (event->type) 
    {
    case Expose :
      XEvent result;
      while (XCheckTypedWindowEvent(me->DisplayId, me->WindowId,
				    Expose, &result))
	{
	// just getting the expose configure event
	event = &result;
	}
      me->GetRenderWindow()->Render();
      break;
      
    case ConfigureNotify : 
      {
      XEvent result;
      while (XCheckTypedWindowEvent(me->DisplayId, me->WindowId,
				    ConfigureNotify, &result))
	{
	// just getting the last configure event
	event = &result;
	}
      if ((((XConfigureEvent *)event)->width != me->Size[0]) ||
	  (((XConfigureEvent *)event)->height != me->Size[1]))
	{
	me->UpdateSize(((XConfigureEvent *)event)->width,
		       ((XConfigureEvent *)event)->height); 

	me->GetRenderWindow()->Render(); 
	}
      }
      break;

    case ButtonPress : 
      {
      me->SetEventPosition(((XButtonEvent*)event)->x,
			   ((XButtonEvent*)event)->y);
      
      switch (((XButtonEvent *)event)->button)
	{
	case Button1 : 
	  if (me->LeftButtonPressMethod) 
	    {
	    (*me->LeftButtonPressMethod)(me->LeftButtonPressMethodArg);
	    }
	  else
	    {
	    me->FindPokedCamera(((XButtonEvent*)event)->x,
				me->Size[1] - ((XButtonEvent*)event)->y);
	    me->StartRotate(); 
	    }
	  break;
	case Button2 : 
	  if (me->MiddleButtonPressMethod) 
	    {
	    (*me->MiddleButtonPressMethod)(me->MiddleButtonPressMethodArg);
	    }
	  else
	    {
	    me->FindPokedCamera(((XButtonEvent*)event)->x,
				me->Size[1] - ((XButtonEvent*)event)->y);
	    me->StartPan(); 
	    }
	  break;
	case Button3 : 
	  if (me->RightButtonPressMethod) 
	    {
	    (*me->RightButtonPressMethod)(me->RightButtonPressMethodArg);
	    }
	  else
	    {
	    me->FindPokedCamera(((XButtonEvent*)event)->x,
				me->Size[1] - ((XButtonEvent*)event)->y);
	    me->StartZoom(); 
	    }
	  break;
	}
      }
      break;

    case ButtonRelease : 
      {
      me->SetEventPosition(((XButtonEvent*)event)->x,
			   ((XButtonEvent*)event)->y);
      switch (((XButtonEvent *)event)->button)
	{
	case Button1 :
	  if (me->LeftButtonReleaseMethod) 
	    {
	    (*me->LeftButtonReleaseMethod)(me->LeftButtonReleaseMethodArg);
	    }
	  else me->EndRotate(); 
	  break;
	case Button2 :
	  if (me->MiddleButtonReleaseMethod) 
	    {
	    (*me->MiddleButtonReleaseMethod)(me->MiddleButtonReleaseMethodArg);
	    }
	  else me->EndPan(); 
	  break;
	case Button3 : 
	  if (me->RightButtonReleaseMethod) 
	    {
	    (*me->RightButtonReleaseMethod)(me->RightButtonReleaseMethodArg);
	    }
	  else me->EndZoom(); 
	  break;
	}
      }
      break;

    case KeyPress :
      {
      KeySym ks;
      static char buffer[20];

      XLookupString((XKeyEvent *)event,buffer,20,&ks,NULL);
      switch (ks)
	{
	case XK_e : 
	  if (me->ExitMethod) (*me->ExitMethod)(me->ExitMethodArg);
	  else Tcl_Exit(1);
	  break;
	case XK_u :
	  if (me->UserMethod) (*me->UserMethod)(me->UserMethodArg);
	  break;
	case XK_r : //reset
	  {
          me->FindPokedRenderer(((XKeyEvent*)event)->x,
			        me->Size[1] - ((XKeyEvent*)event)->y);
	  me->CurrentRenderer->ResetCamera();
	  me->RenderWindow->Render();
          }
	  break;

	case XK_w : //change all actors to wireframe
	  {
	  vtkActorCollection *ac;
	  vtkActor *anActor, *aPart;
	  
          me->FindPokedRenderer(((XKeyEvent*)event)->x,
				me->Size[1] - ((XKeyEvent*)event)->y);
	  ac = me->CurrentRenderer->GetActors();
	  for (ac->InitTraversal(); (anActor = ac->GetNextItem()); )
	    {
            for (anActor->InitPartTraversal();(aPart=anActor->GetNextPart());)
              {
              aPart->GetProperty()->SetRepresentationToWireframe();
              }
	    }
	  
	  me->RenderWindow->Render();
	  }
	  break;

	case XK_s : //change all actors to "surface" or solid
	  {
	  vtkActorCollection *ac;
	  vtkActor *anActor, *aPart;
	  
          me->FindPokedRenderer(((XKeyEvent*)event)->x,
			        me->Size[1] - ((XKeyEvent*)event)->y);
	  ac = me->CurrentRenderer->GetActors();
	  for (ac->InitTraversal(); (anActor = ac->GetNextItem()); )
	    {
            for (anActor->InitPartTraversal();(aPart=anActor->GetNextPart()); )
              {
              aPart->GetProperty()->SetRepresentationToSurface();
              }
	    }
	  
	  me->RenderWindow->Render();
          }
	  break;

	case XK_3 : //3d stereo
	  {
	  // prepare the new window
	  if (me->RenderWindow->GetStereoRender())
	    {
	    if (me->RenderWindow->GetRemapWindow())
	      {
	      me->SetupNewWindow(1);
	      }
	    me->RenderWindow->StereoRenderOff();
	    }
	  else
	    {
	    memcpy(me->PositionBeforeStereo,me->RenderWindow->GetPosition(),
		   sizeof(int)*2);
	    if (me->RenderWindow->GetRemapWindow())
	      {
	      me->SetupNewWindow(1);
	      }
	    me->RenderWindow->StereoRenderOn();
	    }
	  me->RenderWindow->Render();
	  if (me->RenderWindow->GetRemapWindow())
	    {
	    me->FinishSettingUpNewWindow();
	    }
          }
	  break;

	case XK_p : //pick actors
	  {
          me->FindPokedRenderer(((XKeyEvent*)event)->x,
			        me->Size[1] - ((XKeyEvent*)event)->y);
          // Execute start method, if any

          if ( me->StartPickMethod ) 
            (*me->StartPickMethod)(me->StartPickMethodArg);
          me->Picker->Pick(((XButtonEvent*)event)->x,
                             me->Size[1] - ((XButtonEvent*)event)->y, 0.0,
                             me->CurrentRenderer);
          me->HighlightActor(me->Picker->GetAssembly());
          if ( me->EndPickMethod ) 
            (*me->EndPickMethod)(me->EndPickMethodArg);
          }
	  break;
        }
      }
      break;
    }
}

void vtkXRenderWindowInteractorTimer(XtPointer client_data,
				     XtIntervalId *vtkNotUsed(id))
{
  vtkXRenderWindowInteractor *me;
  Window root,child;
  int root_x,root_y;
  int x,y;
  float xf,yf;
  unsigned int keys;

  me = (vtkXRenderWindowInteractor *)client_data;

  if (me->TimerMethod) 
    {
    XQueryPointer(me->DisplayId,me->WindowId,
		  &root,&child,&root_x,&root_y,&x,&y,&keys);
    me->SetEventPosition(x,y);
    (*me->TimerMethod)(me->TimerMethodArg);
    }

  switch (me->State)
    {
    case VTKXI_ROTATE :
      // get the pointer position
      XQueryPointer(me->DisplayId,me->WindowId,
		    &root,&child,&root_x,&root_y,&x,&y,&keys);
      xf = (x - me->Center[0]) * me->DeltaAzimuth;
      yf = ((me->Size[1] - y) - me->Center[1]) * me->DeltaElevation;
      me->CurrentCamera->Azimuth(xf);
      me->CurrentCamera->Elevation(yf);
      me->CurrentCamera->OrthogonalizeViewUp();
      if (me->LightFollowCamera)
	{
	/* get the first light */
	me->CurrentLight->SetPosition(me->CurrentCamera->GetPosition());
	me->CurrentLight->SetFocalPoint(me->CurrentCamera->GetFocalPoint());
	}
      me->RenderWindow->Render();
      Tk_CreateTimerHandler(10,vtkXTclTimerProc,(ClientData)client_data);
      break;
    case VTKXI_PAN :
      {
      float  FPoint[3];
      float *PPoint;
      float  APoint[3];
      float  RPoint[4];

      // get the current focal point and position
      memcpy(FPoint,me->CurrentCamera->GetFocalPoint(),sizeof(float)*3);
      PPoint = me->CurrentCamera->GetPosition();

      // get the pointer position
      XQueryPointer(me->DisplayId,me->WindowId,
		    &root,&child,&root_x,&root_y,&x,&y,&keys);

      APoint[0] = x;
      APoint[1] = me->Size[1] - y;
      APoint[2] = me->FocalDepth;
      me->CurrentRenderer->SetDisplayPoint(APoint);
      me->CurrentRenderer->DisplayToWorld();
      memcpy(RPoint,me->CurrentRenderer->GetWorldPoint(),sizeof(float)*4);
      if (RPoint[3])
	{
	RPoint[0] = RPoint[0]/RPoint[3];
	RPoint[1] = RPoint[1]/RPoint[3];
	RPoint[2] = RPoint[2]/RPoint[3];
	}

      /*
       * Compute a translation vector, moving everything 1/10 
       * the distance to the cursor. (Arbitrary scale factor)
       */
      me->CurrentCamera->SetFocalPoint(
	(FPoint[0]-RPoint[0])/10.0 + FPoint[0],
	(FPoint[1]-RPoint[1])/10.0 + FPoint[1],
	(FPoint[2]-RPoint[2])/10.0 + FPoint[2]);
      me->CurrentCamera->SetPosition(
	(FPoint[0]-RPoint[0])/10.0 + PPoint[0],
	(FPoint[1]-RPoint[1])/10.0 + PPoint[1],
	(FPoint[2]-RPoint[2])/10.0 + PPoint[2]);
      
      me->RenderWindow->Render();
      Tk_CreateTimerHandler(10,vtkXTclTimerProc,(ClientData)client_data);
      }
      break;
    case VTKXI_ZOOM :
      {
      float zoomFactor;
      float *clippingRange;

      // get the pointer position
      XQueryPointer(me->DisplayId,me->WindowId,
		    &root,&child,&root_x,&root_y,&x,&y,&keys);
      yf = ((me->Size[1] - y) - me->Center[1])/(float)me->Center[1];
      zoomFactor = pow(1.1,yf);
      if (me->CurrentCamera->GetParallelProjection())
	{
	me->CurrentCamera->
	  SetParallelScale(me->CurrentCamera->GetParallelScale()/zoomFactor);
	}
      else
	{
	clippingRange = me->CurrentCamera->GetClippingRange();
	me->CurrentCamera->SetClippingRange(clippingRange[0]/zoomFactor,
					    clippingRange[1]/zoomFactor);
	me->CurrentCamera->Dolly(zoomFactor);
	}
      me->RenderWindow->Render();
      Tk_CreateTimerHandler(10,vtkXTclTimerProc,(ClientData)client_data);
      }
      break;
    }
}  


// Description:
// Setup a new window before a WindowRemap
void vtkXRenderWindowInteractor::SetupNewWindow(int Stereo)
{
  vtkXRenderWindow *ren;
  int depth;
  Colormap cmap;
  Visual  *vis;
  int *size;
  int *position;
  int zero_pos[2];
  
  // get the info we need from the RenderingWindow
  ren = (vtkXRenderWindow *)(this->RenderWindow);
  this->DisplayId = ren->GetDisplayId();
  depth   = ren->GetDesiredDepth();
  cmap    = ren->GetDesiredColormap();
  vis     = ren->GetDesiredVisual();
  size    = ren->GetSize();
  position= ren->GetPosition();

  if (Stereo)
    {
    if (this->RenderWindow->GetStereoRender())
      {
      position = this->PositionBeforeStereo;
      }
    else
      {
      zero_pos[0] = 0;
      zero_pos[1] = 0;
      position = zero_pos;
      }
    }
}

// Description:
// Finish setting up a new window after the WindowRemap.
void vtkXRenderWindowInteractor::FinishSettingUpNewWindow()
{
  int *size;

  // free the previous widget
  XSync(this->DisplayId,False);
  this->WindowId = ((vtkXRenderWindow *)(this->RenderWindow))->GetWindowId();
  XSync(this->DisplayId,False);

  XSelectInput(this->DisplayId,this->WindowId,
		    KeyPressMask | ButtonPressMask | ExposureMask |
		    StructureNotifyMask | ButtonReleaseMask);

  size = this->RenderWindow->GetSize();
  this->Size[0] = size[0];
  this->Size[1] = size[1];
}

