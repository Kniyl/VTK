/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkRenderer.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkRenderer - abstract specification for renderers
// .SECTION Description
// vtkRenderer provides an abstract specification for renderers. A renderer
// is an object that controls the rendering process for objects. Rendering
// is the process of converting geometry, a specification for lights, and 
// a camera view into an image. vtkRenderer also performs coordinate 
// transformation between world coordinates, view coordinates (the computer
// graphics rendering coordinate system), and display coordinates (the 
// actual screen coordinates on the display device). Certain advanced 
// rendering features such as two-sided lighting can also be controlled.

// .SECTION See Also
// vtkRenderWindow vtkActor vtkCamera vtkLight vtkVolume

#ifndef __vtkRenderer_h
#define __vtkRenderer_h

#include "vtkViewport.h"

#include "vtkVolumeCollection.h" // Needed for access in inline members
#include "vtkActorCollection.h" // Needed for access in inline members

class vtkRenderWindow;
class vtkVolume;
class vtkCuller;
class vtkActor;
class vtkActor2D;
class vtkCamera;
class vtkLightCollection;
class vtkCullerCollection;
class vtkLight;

class VTK_RENDERING_EXPORT vtkRenderer : public vtkViewport
{
public:
  vtkTypeRevisionMacro(vtkRenderer,vtkViewport);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Create a vtkRenderer with a black background, a white ambient light,
  // two-sided lighting turned on, a viewport of (0,0,1,1), and backface
  // culling turned off.
  static vtkRenderer *New();

  // Description:
  // Add/Remove different types of props to the renderer.
  // These methods are all synonyms to AddProp and RemoveProp.
  // They are here for convenience and backwards compatibility.
  void AddActor(vtkProp *p) {this->AddProp(p);};
  void AddVolume(vtkProp *p) {this->AddProp(p);};
  void RemoveActor(vtkProp *p) {this->Actors->RemoveItem(p);this->RemoveProp(p);};
  void RemoveVolume(vtkProp *p) {this->Volumes->RemoveItem(p);this->RemoveProp(p);};

  // Description:
  // Add a light to the list of lights.
  void AddLight(vtkLight *);

  // Description:
  // Remove a light from the list of lights.
  void RemoveLight(vtkLight *);

  // Description:
  // Return the collection of lights.
  vtkLightCollection *GetLights();
  
  // Description:
  // Create and add a light to renderer.
  void CreateLight(void);
  
  // Description:
  // Create a new Light sutible for use with this type of Renderer.
  // For example, a vtkMesaRenderer should create a vtkMesaLight 
  // in this function.   The default is to just call vtkLight::New.
  virtual vtkLight *MakeLight();

  // Description:
  // Turn on/off two-sided lighting of surfaces. If two-sided lighting is
  // off, then only the side of the surface facing the light(s) will be lit,
  // and the other side dark. If two-sided lighting on, both sides of the 
  // surface will be lit.
  vtkGetMacro(TwoSidedLighting,int);
  vtkSetMacro(TwoSidedLighting,int);
  vtkBooleanMacro(TwoSidedLighting,int);

  // Description:
  // Turn on/off the automatic repositioning of lights as the camera moves.
  // If LightFollowCamera is on, lights that are designated as Headlights
  // or CameraLights will be adjusted to move with this renderer's camera.
  // If LightFollowCamera is off, the lights will not be adjusted.  
  //
  // (Note: In previous versions of vtk, this light-tracking
  // functionality was part of the interactors, not the renderer. For
  // backwards compatibility, the older, more limited interactor
  // behavior is enabled by default. To disable this mode, turn the
  // interactor's LightFollowCamera flag OFF, and leave the renderer's
  // LightFollowCamera flag ON.)
  vtkSetMacro(LightFollowCamera,int);
  vtkGetMacro(LightFollowCamera,int);
  vtkBooleanMacro(LightFollowCamera,int);

  // Description:
  // Turn on/off a flag which disables the automatic light creation capability.
  // Normally in VTK if no lights are associated with the renderer, then a light
  // is automatically created. However, in special circumstances this feature is
  // undesirable, so the following boolean is provided to disable automatic
  // light creation. (Turn AutomaticLightCreation off if you do not want lights
  // to be created.)
  vtkGetMacro(AutomaticLightCreation,int);
  vtkSetMacro(AutomaticLightCreation,int);
  vtkBooleanMacro(AutomaticLightCreation,int);

  // Description:
  // Ask the lights in the scene that are not in world space
  // (for instance, Headlights or CameraLights that are attached to the 
  // camera) to update their geometry to match the active camera.
  virtual int UpdateLightsGeometryToFollowCamera(void);

  // Description:
  // Return the collection of volumes.
  vtkVolumeCollection *GetVolumes();

  // Description:
  // Return any actors in this renderer.
  vtkActorCollection *GetActors();

  // Description:
  // Specify the camera to use for this renderer.
  void SetActiveCamera(vtkCamera *);

  // Description:
  // Get the current camera.
  vtkCamera *GetActiveCamera();

  // Description:
  // Create a new Camera sutible for use with this type of Renderer.
  // For example, a vtkMesaRenderer should create a vtkMesaCamera 
  // in this function.   The default is to just call vtkCamera::New.
  virtual vtkCamera *MakeCamera();

  // Description:
  // Add an culler to the list of cullers.
  void AddCuller(vtkCuller *);

  // Description:
  // Remove an actor from the list of cullers.
  void RemoveCuller(vtkCuller *);

  // Description:
  // Return the collection of cullers.
  vtkCullerCollection *GetCullers();

  // Description:
  // Set the intensity of ambient lighting.
  vtkSetVector3Macro(Ambient,float);
  vtkGetVectorMacro(Ambient,float,3);

  // Description:
  // Set/Get the amount of time this renderer is allowed to spend
  // rendering its scene. This is used by vtkLODActor's.
  vtkSetMacro(AllocatedRenderTime,float);
  virtual float GetAllocatedRenderTime();

  // Description:
  // Get the ratio between allocated time and actual render time.
  // TimeFactor has been taken out of the render process.  
  // It is still computed in case someone finds it useful.
  // It may be taken away in the future.
  virtual float GetTimeFactor();

  // Description:
  // Create an image. This is a superclass method which will in turn 
  // call the DeviceRender method of Subclasses of vtkRenderer
  virtual void Render();

  // Description:
  // Create an image. Subclasses of vtkRenderer must implement this method.
  virtual void DeviceRender() =0;

  // Description:
  // Clear the image to the background color.
  virtual void Clear() {};

  // Description:
  // Returns the number of visible actors.
  int VisibleActorCount();

  // Description:
  // Returns the number of visible volumes.
  int VisibleVolumeCount();

  // Description:
  // Compute the bounding box of all the visible props
  // Used in ResetCamera() and ResetCameraClippingRange() 
  void ComputeVisiblePropBounds( float bounds[6] );

  // Description:
  // Wrapper-friendly version of ComputeVisiblePropBounds 
  float *ComputeVisiblePropBounds();

  // Description:
  // Reset the camera clipping range based on the bounds of the
  // visible actors. This ensures that no props are cut off
  void ResetCameraClippingRange();

  // Description:
  // Reset the camera clipping range based on a bounding box.
  // This method is called from ResetCameraClippingRange()
  void ResetCameraClippingRange( float bounds[6] );
  void ResetCameraClippingRange( float xmin, float xmax, 
                                 float ymin, float ymax, 
                                 float zmin, float zmax);

  // Description:
  // Specify tolerance for near clipping plane distance to the camera
  // as a percentage of the far clipping plane distance.
  vtkSetMacro(NearClippingPlaneTolerance,float);
  vtkGetMacro(NearClippingPlaneTolerance,float);

  // Description:
  // Automatically set up the camera based on the visible actors.
  // The camera will reposition itself to view the center point of the actors,
  // and move along its initial view plane normal (i.e., vector defined from 
  // camera position to focal point) so that all of the actors can be seen.
  void ResetCamera();

  // Description:
  // Automatically set up the camera based on a specified bounding box
  // (xmin,xmax, ymin,ymax, zmin,zmax). Camera will reposition itself so
  // that its focal point is the center of the bounding box, and adjust its
  // distance and position to preserve its initial view plane normal 
  // (i.e., vector defined from camera position to focal point). Note: is 
  // the view plane is parallel to the view up axis, the view up axis will
  // be reset to one of the three coordinate axes.
  void ResetCamera(float bounds[6]);

  // Description:
  // Alternative version of ResetCamera(bounds[6]);
  void ResetCamera(float xmin, float xmax, float ymin, float ymax, 
                   float zmin, float zmax);

  // Description:
  // Specify the rendering window in which to draw. This is automatically set
  // when the renderer is created by MakeRenderer.  The user probably
  // shouldn't ever need to call this method.
  void SetRenderWindow(vtkRenderWindow *);
  vtkRenderWindow *GetRenderWindow() {return this->RenderWindow;};
  virtual vtkWindow *GetVTKWindow();
  
  // Description:
  // Turn on/off using backing store. This may cause the re-rendering
  // time to be slightly slower when the view changes. But it is
  // much faster when the image has not changed, such as during an
  // expose event.
  vtkSetMacro(BackingStore,int);
  vtkGetMacro(BackingStore,int);
  vtkBooleanMacro(BackingStore,int);

  // Description:
  // Turn on/off interactive status.  An interactive renderer is one that 
  // can receive events from an interactor.  Should only be set if
  // there are multiple renderers in the same section of the viewport.
  vtkSetMacro(Interactive,int);
  vtkGetMacro(Interactive,int);
  vtkBooleanMacro(Interactive,int);

  // Description:
  // Set/Get the layer that this renderer belongs to.  This is only used if
  // there are layered renderers.
  vtkSetMacro(Layer, int);
  vtkGetMacro(Layer, int);

  // Description:
  // Returns a boolean indicating if this renderer is transparent.  It is
  // transparent if it is not in the deepest layer of its render window.
  int  Transparent();

  // Description:
  // Convert world point coordinates to view coordinates.
  void WorldToView();

  // Description:
  // Convert view point coordinates to world coordinates.
  void ViewToWorld();
  virtual void ViewToWorld(float &wx, float &wy, float &wz);

  // Description:
  // Convert world point coordinates to view coordinates.
  virtual void WorldToView(float &wx, float &wy, float &wz);

  // Description:
  // Given a pixel location, return the Z value
  float GetZ (int x, int y);

  // Description:
  // Return the MTime of the renderer also considering its ivars.
  unsigned long GetMTime();

  // Description:
  // Get the time required, in seconds, for the last Render call.
  vtkGetMacro( LastRenderTimeInSeconds, float );

  // Description:
  // Should be used internally only during a render
  // Get the number of props that were rendered using a
  // RenderOpaqueGeometry or RenderTranslucentGeometry call.
  // This is used to know if something is in the frame buffer.
  vtkGetMacro( NumberOfPropsRendered, int );

  // Description:
  // Return the prop (via a vtkAssemblyPath) that has the highest z value 
  // at the given x, y position in the viewport.  Basically, the top most 
  // prop that renders the pixel at selectionX, selectionY will be returned.
  // If nothing was picked then NULL is returned.  This method selects from 
  // the renderers Prop list.
  vtkAssemblyPath* PickProp(float selectionX, float selectionY);

protected:
  vtkRenderer();
  ~vtkRenderer();

  // internal method for doing a render for picking purposes
  virtual void PickRender(vtkPropCollection *props);
  virtual void PickGeometry();
  
  vtkCamera *ActiveCamera;
  vtkLight  *CreatedLight;

  vtkLightCollection *Lights;
  vtkCullerCollection *Cullers;

  vtkActorCollection *Actors;
  vtkVolumeCollection *Volumes;
  
  float              Ambient[3];  
  vtkRenderWindow    *RenderWindow;
  float              AllocatedRenderTime;
  float              TimeFactor;
  int                TwoSidedLighting;
  int                AutomaticLightCreation;
  int                BackingStore;
  unsigned char      *BackingImage;
  vtkTimeStamp       RenderTime;

  float              LastRenderTimeInSeconds;

  int                LightFollowCamera;

  // Allocate the time for each prop
  void               AllocateTime();

  // Internal variables indicating the number of props
  // that have been or will be rendered in each category.
  int                NumberOfPropsRendered;

  // A temporary list of props used for culling, and traversal
  // of all props when rendering
  vtkProp            **PropArray;
  int                PropArrayCount;

  // A temporary list used for picking
  vtkAssemblyPath    **PathArray;
  int                PathArrayCount;

  // Indicates if the renderer should receive events from an interactor.
  // Typically only used in conjunction with transparent renderers.
  int                Interactive;

  // Shows what layer this renderer belongs to.  Only of interested when
  // there are layered renderers.
  int                Layer;

  // Holds the result of ComputeVisiblePropBounds so that it is visible from wrapped languages
  float              ComputedVisiblePropBounds[6];

  // Description:
  // Specifies the minimum distance of the near clipping
  // plane as a percentage of the far clipping plane distance.  Values below
  // this threshold are clipped to NearClippingPlaneTolerance*range[1].
  // Note that values which are too small may cause problems on systems
  // with low z-buffer resolution.
  float              NearClippingPlaneTolerance;

  // Description:
  // Ask all props to update and draw any opaque and translucent
  // geometry. This includes both vtkActors and vtkVolumes
  // Returns the number of props that rendered geometry.
  virtual int UpdateGeometry(void);

  // Description:
  // Ask the active camera to do whatever it needs to do prior to rendering.
  // Creates a camera if none found active.
  virtual int UpdateCamera(void);

  // Description:
  // Update the geometry of the lights in the scene that are not in world space
  // (for instance, Headlights or CameraLights that are attached to the camera).
  virtual int UpdateLightGeometry(void);

  // Description:
  // Ask all lights to load themselves into rendering pipeline.
  // This method will return the actual number of lights that were on.
  virtual int UpdateLights(void) {return 0;};
  
private:
  vtkRenderer(const vtkRenderer&);  // Not implemented.
  void operator=(const vtkRenderer&);  // Not implemented.
};

inline vtkLightCollection *vtkRenderer::GetLights() {
  return this->Lights;
}

// Description:
// Get the list of cullers for this renderer.
inline vtkCullerCollection *vtkRenderer::GetCullers(){return this->Cullers;}


#endif
