/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkActor.cxx
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

  Copyright (c) 1993-2002 Ken Martin, Will Schroeder, Bill Lorensen 
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkActor.h"

#include "vtkAssemblyPaths.h"
#include "vtkGraphicsFactory.h"
#include "vtkMapper.h"
#include "vtkPropCollection.h"
#include "vtkProperty.h"
#include "vtkRenderWindow.h"
#include "vtkTexture.h"

#include <math.h>

vtkCxxRevisionMacro(vtkActor, "1.116");

vtkCxxSetObjectMacro(vtkActor,Texture,vtkTexture);


// Creates an actor with the following defaults: origin(0,0,0) 
// position=(0,0,0) scale=(1,1,1) visibility=1 pickable=1 dragable=1
// orientation=(0,0,0). No user defined matrix and no texture map.
vtkActor::vtkActor()
{
  this->Mapper = NULL;
  this->Property = NULL;
  this->BackfaceProperty = NULL;
  this->Texture = NULL;
  
  // The mapper bounds are cache to know when the bounds must be recomputed
  // from the mapper bounds.
  this->MapperBounds[0] = this->MapperBounds[2] = this->MapperBounds[4] = 
  this->MapperBounds[1] = this->MapperBounds[3] = this->MapperBounds[5] = 0.0;

}

vtkActor::~vtkActor()
{
  if ( this->Property != NULL) 
    {
    this->Property->UnRegister(this);
    this->Property = NULL;
    }
  
  if ( this->BackfaceProperty != NULL) 
    {
    this->BackfaceProperty->UnRegister(this);
    this->BackfaceProperty = NULL;
    }

  if (this->Mapper)
    {
    this->Mapper->UnRegister(this);
    this->Mapper = NULL;
    }
  this->SetTexture(NULL);
}

// Shallow copy of an actor.
void vtkActor::ShallowCopy(vtkProp *prop)
{
  vtkActor *a = vtkActor::SafeDownCast(prop);
  if ( a != NULL )
    {
    this->SetMapper(a->GetMapper());
    this->SetProperty(a->GetProperty());
    this->SetBackfaceProperty(a->GetBackfaceProperty());
    this->SetTexture(a->GetTexture());
    }

  // Now do superclass
  this->vtkProp3D::ShallowCopy(prop);
}

// return the correct type of Actor 
vtkActor *vtkActor::New()
{
  // First try to create the object from the vtkGraphicsFactory
  vtkObject* ret = vtkGraphicsFactory::CreateInstance("vtkActor");
  return (vtkActor*)ret;
}

void vtkActor::GetActors(vtkPropCollection *ac)
{
  ac->AddItem(this);
}

// should be called from the render methods only
int vtkActor::GetIsOpaque()
{
  if (this->Property->GetOpacity() >= 1.0)
    {
    if (this->Texture && this->Texture->GetInput())
      {
      this->Texture->GetInput()->UpdateInformation();
      this->Texture->GetInput()->SetUpdateExtent(this->Texture->GetInput()->GetWholeExtent());
      this->Texture->GetInput()->PropagateUpdateExtent();
      this->Texture->GetInput()->TriggerAsynchronousUpdate();
      this->Texture->GetInput()->UpdateData();
      if (this->Texture->GetInput()->GetPointData()->GetScalars() == NULL)
        { // Handle gracefully. What should it return? 
        return 1;
        }
      if (this->Texture->GetInput()->GetPointData()->GetScalars()->GetNumberOfComponents()%2)
        {
        return 1;
        }
      }
    else
      {
      return 1;
      }
    }
  return 0;
}


// This causes the actor to be rendered. It in turn will render the actor's
// property, texture map and then mapper. If a property hasn't been 
// assigned, then the actor will create one automatically. Note that a 
// side effect of this method is that the visualization network is updated.
int vtkActor::RenderOpaqueGeometry(vtkViewport *vp)
{
  int          renderedSomething = 0; 
  vtkRenderer  *ren = (vtkRenderer *)vp;

  if ( ! this->Mapper )
    {
    return 0;
    }

  // make sure we have a property
  if (!this->Property)
    {
    // force creation of a property
    this->GetProperty();
    }

  // is this actor opaque ?
  if (this->GetIsOpaque())
    {
    this->Property->Render(this, ren);

    // render the backface property
    if (this->BackfaceProperty)
      {
      this->BackfaceProperty->BackfaceRender(this, ren);
      }
    
    // render the texture 
    if (this->Texture)
      {
      this->Texture->Render(ren);
      }
    this->Render(ren,this->Mapper);
    this->EstimatedRenderTime += this->Mapper->GetTimeToDraw();

    renderedSomething = 1;
    }

  return renderedSomething;
}

int vtkActor::RenderTranslucentGeometry(vtkViewport *vp)
{
  int          renderedSomething = 0; 
  vtkRenderer *ren = (vtkRenderer *)vp;

  if ( ! this->Mapper )
    {
    return 0;
    }

  // make sure we have a property
  if (!this->Property)
    {
    // force creation of a property
    this->GetProperty();
    }

  // is this actor opaque ?
  if (!this->GetIsOpaque())
    {
    this->Property->Render(this, ren);

    // render the backface property
    if (this->BackfaceProperty)
      {
      this->BackfaceProperty->BackfaceRender(this, ren);
      }
    
    // render the texture 
    if (this->Texture)
      {
      this->Texture->Render(ren);
      }
    this->Render(ren,this->Mapper);
    this->EstimatedRenderTime += this->Mapper->GetTimeToDraw();

    renderedSomething = 1;
    }

  return renderedSomething;
}

void vtkActor::ReleaseGraphicsResources(vtkWindow *win)
{
  vtkRenderWindow *renWin = (vtkRenderWindow *)win;

  // pass this information onto the mapper
  if (this->Mapper)
    {
    this->Mapper->ReleaseGraphicsResources(renWin);
    }

  // pass this information onto the texture
  if (this->Texture)
    {
    this->Texture->ReleaseGraphicsResources(renWin);
    }
}

void vtkActor::SetProperty(vtkProperty *lut)
{
  if ( this->Property == lut) 
    {
    return;
    }
  if ( this->Property != NULL) 
    {
    this->Property->UnRegister(this);
    this->Property = NULL;
    }
  if ( lut != NULL) 
    {
    lut->Register(this);
    }
  
  this->Property = lut;
  this->Modified();
}

vtkProperty* vtkActor::MakeProperty()
{
  return vtkProperty::New();
}


vtkProperty *vtkActor::GetProperty()
{
  if ( this->Property == NULL )
    {
    vtkProperty *p = this->MakeProperty();
    this->SetProperty(p);
    p->Delete();
    }
  return this->Property;
}

void vtkActor::SetBackfaceProperty(vtkProperty *lut)
{
  if ( this->BackfaceProperty == lut) 
    {
    return;
    }
  if ( this->BackfaceProperty != NULL) 
    {
    this->BackfaceProperty->UnRegister(this);
    this->BackfaceProperty = NULL;
    }
  if ( lut != NULL) 
    {
    lut->Register(this);
    }
  
  this->BackfaceProperty = lut;
  this->Modified();
}

vtkProperty *vtkActor::GetBackfaceProperty()
{
  return this->BackfaceProperty;
}

// Get the bounds for this Actor as (Xmin,Xmax,Ymin,Ymax,Zmin,Zmax).
float *vtkActor::GetBounds()
{
  int i,n;
  float *bounds, bbox[24], *fptr;

  vtkDebugMacro( << "Getting Bounds" );

  // get the bounds of the Mapper if we have one
  if (!this->Mapper)
    {
    return this->Bounds;
    }

  bounds = this->Mapper->GetBounds();
  // Check for the special case when the mapper's bounds are unknown
  if (!bounds)
    {
    return bounds;
    }

  // Check for the special case when the actor is empty.
  if (bounds[0] > bounds[1])
    { 
    memcpy( this->MapperBounds, bounds, 6*sizeof(float) );
    this->Bounds[0] = this->Bounds[2] = this->Bounds[4] = VTK_LARGE_FLOAT;
    this->Bounds[1] = this->Bounds[3] = this->Bounds[5] = -VTK_LARGE_FLOAT;
    this->BoundsMTime.Modified();
    return this->Bounds;
    }

  // Check if we have cached values for these bounds - we cache the
  // values returned by this->Mapper->GetBounds() and we store the time
  // of caching. If the values returned this time are different, or
  // the modified time of this class is newer than the cached time,
  // then we need to rebuild.
  if ( ( memcmp( this->MapperBounds, bounds, 6*sizeof(float) ) != 0 ) ||
       ( this->GetMTime() > this->BoundsMTime ) )
    {
    vtkDebugMacro( << "Recomputing bounds..." );

    memcpy( this->MapperBounds, bounds, 6*sizeof(float) );

    // fill out vertices of a bounding box
    bbox[ 0] = bounds[1]; bbox[ 1] = bounds[3]; bbox[ 2] = bounds[5];
    bbox[ 3] = bounds[1]; bbox[ 4] = bounds[2]; bbox[ 5] = bounds[5];
    bbox[ 6] = bounds[0]; bbox[ 7] = bounds[2]; bbox[ 8] = bounds[5];
    bbox[ 9] = bounds[0]; bbox[10] = bounds[3]; bbox[11] = bounds[5];
    bbox[12] = bounds[1]; bbox[13] = bounds[3]; bbox[14] = bounds[4];
    bbox[15] = bounds[1]; bbox[16] = bounds[2]; bbox[17] = bounds[4];
    bbox[18] = bounds[0]; bbox[19] = bounds[2]; bbox[20] = bounds[4];
    bbox[21] = bounds[0]; bbox[22] = bounds[3]; bbox[23] = bounds[4];
  
    // save the old transform
    this->Transform->Push(); 
    this->Transform->SetMatrix(this->GetMatrix());

    // and transform into actors coordinates
    fptr = bbox;
    for (n = 0; n < 8; n++) 
      {
      this->Transform->TransformPoint(fptr,fptr);
      fptr += 3;
      }
  
    this->Transform->Pop();  
  
    // now calc the new bounds
    this->Bounds[0] = this->Bounds[2] = this->Bounds[4] = VTK_LARGE_FLOAT;
    this->Bounds[1] = this->Bounds[3] = this->Bounds[5] = -VTK_LARGE_FLOAT;
    for (i = 0; i < 8; i++)
      {
      for (n = 0; n < 3; n++)
        {
        if (bbox[i*3+n] < this->Bounds[n*2]) 
          {
          this->Bounds[n*2] = bbox[i*3+n];
          }
        if (bbox[i*3+n] > this->Bounds[n*2+1]) 
          {
          this->Bounds[n*2+1] = bbox[i*3+n];
          }
        }
      }
    this->BoundsMTime.Modified();
    }

  return this->Bounds;
}

unsigned long int vtkActor::GetMTime()
{
  unsigned long mTime=this->vtkObject::GetMTime();
  unsigned long time;

  if ( this->Property != NULL )
    {
    time = this->Property->GetMTime();
    mTime = ( time > mTime ? time : mTime );
    }

  if ( this->BackfaceProperty != NULL )
    {
    time = this->BackfaceProperty->GetMTime();
    mTime = ( time > mTime ? time : mTime );
    }

  if ( this->UserMatrix != NULL )
    {
    time = this->UserMatrix->GetMTime();
    mTime = ( time > mTime ? time : mTime );
    }

  if ( this->UserTransform != NULL )
    {
    time = this->UserTransform->GetMTime();
    mTime = ( time > mTime ? time : mTime );
    }

  if ( this->Texture != NULL )
    {
    time = this->Texture->GetMTime();
    mTime = ( time > mTime ? time : mTime );
    }

  return mTime;
}

unsigned long int vtkActor::GetRedrawMTime()
{
  unsigned long mTime=this->GetMTime();
  unsigned long time;

  if ( this->Mapper != NULL )
    {
    time = this->Mapper->GetMTime();
    mTime = ( time > mTime ? time : mTime );
    if (this->GetMapper()->GetInput() != NULL)
      {
      this->GetMapper()->GetInput()->Update();
      time = this->Mapper->GetInput()->GetMTime();
      mTime = ( time > mTime ? time : mTime );
      }
    }

  return mTime;
}

void vtkActor::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  if ( this->Mapper )
    {
    os << indent << "Mapper:\n";
    this->Mapper->PrintSelf(os,indent.GetNextIndent());
    }
  else
    {
    os << indent << "Mapper: (none)\n";
    }

  if ( this->Property )
    {
    os << indent << "Property:\n";
    this->Property->PrintSelf(os,indent.GetNextIndent());
    }
  else
    {
    os << indent << "Property: (none)\n";
    }

  if ( this->Texture )
    {
    os << indent << "Texture: " << this->Texture << "\n";
    }
  else
    {
    os << indent << "Texture: (none)\n";
    }

}

// Compatibility methods...to be deprecated in the future.
#include "vtkAssemblyNode.h"
void vtkActor::InitPartTraversal()
{
  this->InitPathTraversal();
}

vtkActor *vtkActor::GetNextPart()
{
  vtkAssemblyPath *path = this->GetNextPath();
  if ( !path )
    {
    return NULL;
    }
  else
    {
    vtkAssemblyNode *node = path->GetLastNode();
    if ( node && node->GetProp()->IsA("vtkActor") )
      {
      return (vtkActor *)node->GetProp();
      }
    }
  return NULL;
}

int vtkActor::GetNumberOfParts()
{
  return this->GetNumberOfPaths();
}


void vtkActor::SetMapper(vtkMapper *args)
{
  vtkDebugMacro(<< this->GetClassName() << " (" << this         
  << "): setting Mapper to " << args );     
  if (this->Mapper != args)                                       
    {                                                           
    if (this->Mapper != NULL) 
      { 
      this->Mapper->RemoveConsumer(this);
      this->Mapper->UnRegister(this); 
      }   
    this->Mapper = args;                                          
    if (this->Mapper != NULL) 
      { 
      this->Mapper->Register(this); 
      this->Mapper->AddConsumer(this);
      }     
    this->Modified();                                           
    }                                                           
}
