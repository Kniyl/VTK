/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkRenderWindowInteractor.cxx
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
#ifdef _WIN32
#include "vtkWin32RenderWindowInteractor.h"
#else
#include "vtkXRenderWindowInteractor.h"
#endif
#include "vtkActor.h"
#include "vtkCellPicker.h"

// Construct object so that light follows camera motion.
vtkRenderWindowInteractor::vtkRenderWindowInteractor()
{
  this->RenderWindow    = NULL;
  this->CurrentCamera   = NULL;
  this->CurrentLight    = NULL;
  this->CurrentRenderer = NULL;

  this->LightFollowCamera = 1;
  this->Initialized = 0;
  this->DesiredUpdateRate = 15;
  // default limit is 3 hours per frame
  this->StillUpdateRate = 0.0001;
  
  this->SelfCreatedPicker = 0;
  this->Picker = this->CreateDefaultPicker();
  this->OutlineActor = NULL;
  this->OutlineMapper = vtkPolyDataMapper::New();  
  this->OutlineMapper->SetInput(this->Outline.GetOutput());
  this->PickedRenderer = NULL;
  this->CurrentActor = NULL;

  this->StartPickMethod = NULL;
  this->StartPickMethodArgDelete = NULL;
  this->StartPickMethodArg = NULL;
  this->EndPickMethod = NULL;
  this->EndPickMethodArgDelete = NULL;
  this->EndPickMethodArg = NULL;
  this->UserMethod = NULL;
  this->UserMethodArgDelete = NULL;
  this->UserMethodArg = NULL;
  this->ExitMethod = NULL;
  this->ExitMethodArgDelete = NULL;
  this->ExitMethodArg = NULL;

  this->TimerMethod = NULL;
  this->TimerMethodArgDelete = NULL;
  this->TimerMethodArg = NULL;

  this->LeftButtonPressMethod = NULL;
  this->LeftButtonPressMethodArgDelete = NULL;
  this->LeftButtonPressMethodArg = NULL;
  this->LeftButtonReleaseMethod = NULL;
  this->LeftButtonReleaseMethodArgDelete = NULL;
  this->LeftButtonReleaseMethodArg = NULL;

  this->MiddleButtonPressMethod = NULL;
  this->MiddleButtonPressMethodArgDelete = NULL;
  this->MiddleButtonPressMethodArg = NULL;
  this->MiddleButtonReleaseMethod = NULL;
  this->MiddleButtonReleaseMethodArgDelete = NULL;
  this->MiddleButtonReleaseMethodArg = NULL;

  this->RightButtonPressMethod = NULL;
  this->RightButtonPressMethodArgDelete = NULL;
  this->RightButtonPressMethodArg = NULL;
  this->RightButtonReleaseMethod = NULL;
  this->RightButtonReleaseMethodArgDelete = NULL;
  this->RightButtonReleaseMethodArg = NULL;

  this->EventPosition[0] = 0;
  this->EventPosition[1] = 0;
}

vtkRenderWindowInteractor::~vtkRenderWindowInteractor()
{
  if ( this->OutlineActor ) this->OutlineActor->Delete();
  if ( this->OutlineMapper ) this->OutlineMapper->Delete();
  if ( this->SelfCreatedPicker && this->Picker) this->Picker->Delete();

  // delete the current arg if there is one and a delete meth
  if ((this->UserMethodArg)&&(this->UserMethodArgDelete))
    {
    (*this->UserMethodArgDelete)(this->UserMethodArg);
    }
  if ((this->ExitMethodArg)&&(this->ExitMethodArgDelete))
    {
    (*this->ExitMethodArgDelete)(this->ExitMethodArg);
    }
  if ((this->StartPickMethodArg)&&(this->StartPickMethodArgDelete))
    {
    (*this->StartPickMethodArgDelete)(this->StartPickMethodArg);
    }
  if ((this->EndPickMethodArg)&&(this->EndPickMethodArgDelete))
    {
    (*this->EndPickMethodArgDelete)(this->EndPickMethodArg);
    }
  if ((this->TimerMethodArg)&&(this->TimerMethodArgDelete))
    {
    (*this->TimerMethodArgDelete)(this->TimerMethodArg);
    }
  if ((this->LeftButtonPressMethodArg)&&(this->LeftButtonPressMethodArgDelete))
    {
    (*this->LeftButtonPressMethodArgDelete)(this->LeftButtonPressMethodArg);
    }
  if ((this->LeftButtonReleaseMethodArg)&&
      (this->LeftButtonReleaseMethodArgDelete))
    {
    (*this->LeftButtonReleaseMethodArgDelete)
      (this->LeftButtonReleaseMethodArg);
    }
  if ((this->MiddleButtonPressMethodArg)&&
      (this->MiddleButtonPressMethodArgDelete))
    {
    (*this->MiddleButtonPressMethodArgDelete)
      (this->MiddleButtonPressMethodArg);
    }
  if ((this->MiddleButtonReleaseMethodArg)&&
      (this->MiddleButtonReleaseMethodArgDelete))
    {
    (*this->MiddleButtonReleaseMethodArgDelete)
      (this->MiddleButtonReleaseMethodArg);
    }
  if ((this->RightButtonPressMethodArg)&&
      (this->RightButtonPressMethodArgDelete))
    {
    (*this->RightButtonPressMethodArgDelete)(this->RightButtonPressMethodArg);
    }
  if ((this->RightButtonReleaseMethodArg)&&
      (this->RightButtonReleaseMethodArgDelete))
    {
    (*this->RightButtonReleaseMethodArgDelete)
      (this->RightButtonReleaseMethodArg);
    }
}

vtkRenderWindowInteractor *vtkRenderWindowInteractor::New()
{
#ifdef _WIN32
  return vtkWin32RenderWindowInteractor::New();
#else
  return vtkXRenderWindowInteractor::New();
#endif  
}

void vtkRenderWindowInteractor::SetRenderWindow(vtkRenderWindow *aren)
{
  this->RenderWindow = aren;
  if (this->RenderWindow->GetInteractor() != this)
    {
    this->RenderWindow->SetInteractor(this);
    }
}
void vtkRenderWindowInteractor::FindPokedRenderer(int x,int y)
{
  vtkRendererCollection *rc;
  vtkRenderer *aren;

  this->CurrentRenderer = NULL;

  rc = this->RenderWindow->GetRenderers();
  
  for (rc->InitTraversal(); 
       ((aren = rc->GetNextItem())&&(!this->CurrentRenderer));)
    {
    if (aren->IsInViewport(x,y))
      {
      this->CurrentRenderer = aren;
      }
    }
  
  // we must have a value 
  if (this->CurrentRenderer == NULL)
    {
    rc->InitTraversal();
    aren = rc->GetNextItem();
    this->CurrentRenderer = aren;
    }
}

void  vtkRenderWindowInteractor::FindPokedCamera(int x,int y)
{
  float *vp;
  vtkLightCollection *lc;

  this->FindPokedRenderer(x,y);
  vp = this->CurrentRenderer->GetViewport();

  this->CurrentCamera = this->CurrentRenderer->GetActiveCamera();  
  memcpy(this->Center,this->CurrentRenderer->GetCenter(),sizeof(int)*2);
  this->DeltaElevation = -20.0/((vp[3] - vp[1])*this->Size[1]);
  this->DeltaAzimuth = -20.0/((vp[2] - vp[0])*this->Size[0]);

  // as a side effect also set the light 
  // in case they are using light follow camera 
  lc = this->CurrentRenderer->GetLights();
  lc->InitTraversal();
  this->CurrentLight = lc->GetNextItem();
}

// When pick action successfully selects actor, this method highlights the 
// actor appropriately. Currently this is done by placing a bounding box
// around the actor.
void vtkRenderWindowInteractor::HighlightActor(vtkActor *actor)
{
  if ( ! this->OutlineActor )
    {
    // have to defer creation to get right type
    this->OutlineActor = vtkActor::New();
    this->OutlineActor->PickableOff();
    this->OutlineActor->DragableOff();
    this->OutlineActor->SetMapper(this->OutlineMapper);
    this->OutlineActor->GetProperty()->SetColor(1.0,1.0,1.0);
    this->OutlineActor->GetProperty()->SetAmbient(1.0);
    this->OutlineActor->GetProperty()->SetDiffuse(0.0);
    }

  if ( this->PickedRenderer ) 
    this->PickedRenderer->RemoveActor(OutlineActor);

  if ( ! actor )
    {
    this->PickedRenderer = NULL;
    }
  else 
    {
    this->PickedRenderer = this->CurrentRenderer;
    this->CurrentRenderer->AddActor(OutlineActor);
    this->Outline.SetBounds(actor->GetBounds());
    this->CurrentActor = actor;
    }
  this->RenderWindow->Render();
}

// Specify a method to be executed prior to the pick operation.
void vtkRenderWindowInteractor::SetStartPickMethod(void (*f)(void *), void *arg)
{
  if ( f != this->StartPickMethod || arg != this->StartPickMethodArg )
    {
    // delete the current arg if there is one and a delete meth
    if ((this->StartPickMethodArg)&&(this->StartPickMethodArgDelete))
      {
      (*this->StartPickMethodArgDelete)(this->StartPickMethodArg);
      }
    this->StartPickMethod = f;
    this->StartPickMethodArg = arg;
    this->Modified();
    }
}

// Specify a method to be executed after the pick operation.
void vtkRenderWindowInteractor::SetEndPickMethod(void (*f)(void *), void *arg)
{
  if ( f != this->EndPickMethod || arg != this->EndPickMethodArg )
    {
    // delete the current arg if there is one and a delete meth
    if ((this->EndPickMethodArg)&&(this->EndPickMethodArgDelete))
      {
      (*this->EndPickMethodArgDelete)(this->EndPickMethodArg);
      }
    this->EndPickMethod = f;
    this->EndPickMethodArg = arg;
    this->Modified();
    }
}

// Set the object used to perform pick operations. You can use this to 
// control what type of data is picked.
void vtkRenderWindowInteractor::SetPicker(vtkPicker *picker)
{
  if ( this->Picker != picker ) 
    {
    if ( this->SelfCreatedPicker ) this->Picker->Delete();
    this->SelfCreatedPicker = 0;
    this->Picker = picker;
    this->Modified();
    }
}

vtkPicker *vtkRenderWindowInteractor::CreateDefaultPicker()
{
  if ( this->SelfCreatedPicker ) this->Picker->Delete();
  this->SelfCreatedPicker = 1;
  return vtkCellPicker::New();
}

// Set the user method. This method is invoked on a <u> keypress.
void vtkRenderWindowInteractor::SetUserMethod(void (*f)(void *), void *arg)
{
  if ( f != this->UserMethod || arg != this->UserMethodArg )
    {
    // delete the current arg if there is one and a delete meth
    if ((this->UserMethodArg)&&(this->UserMethodArgDelete))
      {
      (*this->UserMethodArgDelete)(this->UserMethodArg);
      }
    this->UserMethod = f;
    this->UserMethodArg = arg;
    this->Modified();
    }
}

// Called when a void* argument is being discarded.  Lets the user free it.
void vtkRenderWindowInteractor::SetUserMethodArgDelete(void (*f)(void *))
{
  if ( f != this->UserMethodArgDelete)
    {
    this->UserMethodArgDelete = f;
    this->Modified();
    }
}

// Set the exit method. This method is invoked on a <e> keypress.
void vtkRenderWindowInteractor::SetExitMethod(void (*f)(void *), void *arg)
{
  if ( f != this->ExitMethod || arg != this->ExitMethodArg )
    {
    // delete the current arg if there is one and a delete meth
    if ((this->ExitMethodArg)&&(this->ExitMethodArgDelete))
      {
      (*this->ExitMethodArgDelete)(this->ExitMethodArg);
      }
    this->ExitMethod = f;
    this->ExitMethodArg = arg;
    this->Modified();
    }
}

// Called when a void* argument is being discarded.  Lets the user free it.
void vtkRenderWindowInteractor::SetExitMethodArgDelete(void (*f)(void *))
{
  if ( f != this->ExitMethodArgDelete)
    {
    this->ExitMethodArgDelete = f;
    this->Modified();
    }
}

// Set the exit method. This method is invoked during rotate/zoom/pan
void vtkRenderWindowInteractor::SetTimerMethod(void (*f)(void *), void *arg)
{
  if ( f != this->TimerMethod || arg != this->TimerMethodArg )
    {
    // delete the current arg if there is one and a delete meth
    if ((this->TimerMethodArg)&&(this->TimerMethodArgDelete))
      {
      (*this->TimerMethodArgDelete)(this->TimerMethodArg);
      }
    this->TimerMethod = f;
    this->TimerMethodArg = arg;
    this->Modified();
    }
}

// Called when a void* argument is being discarded.  Lets the user free it.
void vtkRenderWindowInteractor::SetTimerMethodArgDelete(void (*f)(void *))
{
  if ( f != this->TimerMethodArgDelete)
    {
    this->TimerMethodArgDelete = f;
    this->Modified();
    }
}

// Set the exit method. This method is invoked on a <e> keypress.
void vtkRenderWindowInteractor::SetLeftButtonPressMethod(void (*f)(void *), void *arg)
{
  if ( f != this->LeftButtonPressMethod || arg != this->LeftButtonPressMethodArg )
    {
    // delete the current arg if there is one and a delete meth
    if ((this->LeftButtonPressMethodArg)&&(this->LeftButtonPressMethodArgDelete))
      {
      (*this->LeftButtonPressMethodArgDelete)(this->LeftButtonPressMethodArg);
      }
    this->LeftButtonPressMethod = f;
    this->LeftButtonPressMethodArg = arg;
    this->Modified();
    }
}

// Called when a void* argument is being discarded.  Lets the user free it.
void vtkRenderWindowInteractor::SetLeftButtonPressMethodArgDelete(void (*f)(void *))
{
  if ( f != this->LeftButtonPressMethodArgDelete)
    {
    this->LeftButtonPressMethodArgDelete = f;
    this->Modified();
    }
}

// Set the exit method. This method is invoked on a <e> keyrelease.
void vtkRenderWindowInteractor::SetLeftButtonReleaseMethod(void (*f)(void *), void *arg)
{
  if ( f != this->LeftButtonReleaseMethod || arg != this->LeftButtonReleaseMethodArg )
    {
    // delete the current arg if there is one and a delete meth
    if ((this->LeftButtonReleaseMethodArg)&&(this->LeftButtonReleaseMethodArgDelete))
      {
      (*this->LeftButtonReleaseMethodArgDelete)(this->LeftButtonReleaseMethodArg);
      }
    this->LeftButtonReleaseMethod = f;
    this->LeftButtonReleaseMethodArg = arg;
    this->Modified();
    }
}

// Called when a void* argument is being discarded.  Lets the user free it.
void vtkRenderWindowInteractor::SetLeftButtonReleaseMethodArgDelete(void (*f)(void *))
{
  if ( f != this->LeftButtonReleaseMethodArgDelete)
    {
    this->LeftButtonReleaseMethodArgDelete = f;
    this->Modified();
    }
}

// Set the exit method. This method is invoked on a <e> keypress.
void vtkRenderWindowInteractor::SetMiddleButtonPressMethod(void (*f)(void *), void *arg)
{
  if ( f != this->MiddleButtonPressMethod || arg != this->MiddleButtonPressMethodArg )
    {
    // delete the current arg if there is one and a delete meth
    if ((this->MiddleButtonPressMethodArg)&&(this->MiddleButtonPressMethodArgDelete))
      {
      (*this->MiddleButtonPressMethodArgDelete)(this->MiddleButtonPressMethodArg);
      }
    this->MiddleButtonPressMethod = f;
    this->MiddleButtonPressMethodArg = arg;
    this->Modified();
    }
}

// Called when a void* argument is being discarded.  Lets the user free it.
void vtkRenderWindowInteractor::SetMiddleButtonPressMethodArgDelete(void (*f)(void *))
{
  if ( f != this->MiddleButtonPressMethodArgDelete)
    {
    this->MiddleButtonPressMethodArgDelete = f;
    this->Modified();
    }
}

// Set the exit method. This method is invoked on a <e> keyrelease.
void vtkRenderWindowInteractor::SetMiddleButtonReleaseMethod(void (*f)(void *), void *arg)
{
  if ( f != this->MiddleButtonReleaseMethod || arg != this->MiddleButtonReleaseMethodArg )
    {
    // delete the current arg if there is one and a delete meth
    if ((this->MiddleButtonReleaseMethodArg)&&(this->MiddleButtonReleaseMethodArgDelete))
      {
      (*this->MiddleButtonReleaseMethodArgDelete)(this->MiddleButtonReleaseMethodArg);
      }
    this->MiddleButtonReleaseMethod = f;
    this->MiddleButtonReleaseMethodArg = arg;
    this->Modified();
    }
}

// Called when a void* argument is being discarded.  Lets the user free it.
void vtkRenderWindowInteractor::SetMiddleButtonReleaseMethodArgDelete(void (*f)(void *))
{
  if ( f != this->MiddleButtonReleaseMethodArgDelete)
    {
    this->MiddleButtonReleaseMethodArgDelete = f;
    this->Modified();
    }
}

// Set the exit method. This method is invoked on a <e> keypress.
void vtkRenderWindowInteractor::SetRightButtonPressMethod(void (*f)(void *), void *arg)
{
  if ( f != this->RightButtonPressMethod || arg != this->RightButtonPressMethodArg )
    {
    // delete the current arg if there is one and a delete meth
    if ((this->RightButtonPressMethodArg)&&(this->RightButtonPressMethodArgDelete))
      {
      (*this->RightButtonPressMethodArgDelete)(this->RightButtonPressMethodArg);
      }
    this->RightButtonPressMethod = f;
    this->RightButtonPressMethodArg = arg;
    this->Modified();
    }
}

// Called when a void* argument is being discarded.  Lets the user free it.
void vtkRenderWindowInteractor::SetRightButtonPressMethodArgDelete(void (*f)(void *))
{
  if ( f != this->RightButtonPressMethodArgDelete)
    {
    this->RightButtonPressMethodArgDelete = f;
    this->Modified();
    }
}

// Set the exit method. This method is invoked on a <e> keyrelease.
void vtkRenderWindowInteractor::SetRightButtonReleaseMethod(void (*f)(void *), void *arg)
{
  if ( f != this->RightButtonReleaseMethod || arg != this->RightButtonReleaseMethodArg )
    {
    // delete the current arg if there is one and a delete meth
    if ((this->RightButtonReleaseMethodArg)&&(this->RightButtonReleaseMethodArgDelete))
      {
      (*this->RightButtonReleaseMethodArgDelete)(this->RightButtonReleaseMethodArg);
      }
    this->RightButtonReleaseMethod = f;
    this->RightButtonReleaseMethodArg = arg;
    this->Modified();
    }
}

// Called when a void* argument is being discarded.  Lets the user free it.
void vtkRenderWindowInteractor::SetRightButtonReleaseMethodArgDelete(void (*f)(void *))
{
  if ( f != this->RightButtonReleaseMethodArgDelete)
    {
    this->RightButtonReleaseMethodArgDelete = f;
    this->Modified();
    }
}

// Called when a void* argument is being discarded.  Lets the user free it.
void vtkRenderWindowInteractor::SetStartPickMethodArgDelete(void (*f)(void *))
{
  if ( f != this->StartPickMethodArgDelete)
    {
    this->StartPickMethodArgDelete = f;
    this->Modified();
    }
}
// Called when a void* argument is being discarded.  Lets the user free it.
void vtkRenderWindowInteractor::SetEndPickMethodArgDelete(void (*f)(void *))
{
  if ( f != this->EndPickMethodArgDelete)
    {
    this->EndPickMethodArgDelete = f;
    this->Modified();
    }
}

void vtkRenderWindowInteractor::PrintSelf(ostream& os, vtkIndent indent)
{
  vtkObject::PrintSelf(os,indent);

  os << indent << "RenderWindow:    " << this->RenderWindow << "\n";
  os << indent << "CurrentCamera:   " << this->CurrentCamera << "\n";
  os << indent << "CurrentLight:    " << this->CurrentLight << "\n";
  os << indent << "CurrentRenderer: " << this->CurrentRenderer << "\n";
  if ( this->Picker )
    {
    os << indent << "Picker: " << this->Picker << "\n";
    }
  else
    {
    os << indent << "Picker: (none)\n";
    }
  os << indent << "LightFollowCamera: " << (this->LightFollowCamera ? "On\n" : "Off\n");
  os << indent << "DesiredUpdateRate: " << this->DesiredUpdateRate << "\n";
  os << indent << "StillUpdateRate: " << this->StillUpdateRate << "\n";
  os << indent << "Initialized: " << this->Initialized << "\n";
  os << indent << "EventPosition: " << "( " << this->EventPosition[0] <<
     ", " << this->EventPosition[1] << " )\n";
}

