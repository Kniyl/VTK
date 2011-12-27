/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkBrokenLineWidget.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkBrokenLineWidget - 3D widget for manipulating a broken line
// .SECTION Description
// This 3D widget defines a broken line that can be interactively placed in a
// scene. The broken line has handles, the number of which can be changed, plus it
// can be picked on the broken line itself to translate or rotate it in the scene.
// A nice feature of the object is that the vtkBrokenLineWidget, like any 3D
// widget, will work with the current interactor style. That is, if
// vtkBrokenLineWidget does not handle an event, then all other registered
// observers (including the interactor style) have an opportunity to process
// the event. Otherwise, the vtkBrokenLineWidget will terminate the processing of
// the event that it handles.
//
// To use this object, just invoke SetInteractor() with the argument of the
// method a vtkRenderWindowInteractor.  You may also wish to invoke
// "PlaceWidget()" to initially position the widget. The interactor will act
// normally until the "i" key (for "interactor") is pressed, at which point the
// vtkBrokenLineWidget will appear. (See superclass documentation for information
// about changing this behavior.) Events that occur outside of the widget
// (i.e., no part of the widget is picked) are propagated to any other
// registered obsevers (such as the interaction style).  Turn off the widget
// by pressing the "i" key again (or invoke the Off() method).
//
// The button actions and key modifiers are as follows for controlling the
// widget:
// 1) left button down on and drag one of the spherical handles to change the
// shape of the broken line: the handles act as "control points".
// 2) left button or middle button down on a line segment forming the broken line
// allows uniform translation of the widget.
// 3) ctrl + middle button down on the widget enables spinning of the widget
// about its center.
// 4) right button down on the widget enables scaling of the widget. By moving
// the mouse "up" the render window the broken line will be made bigger; by moving
// "down" the render window the widget will be made smaller.
// 5) ctrl key + right button down on any handle will erase it providing there
// will be two or more points remaining to form a broken line.
// 6) shift key + right button down on any line segment will insert a handle
// onto the broken line at the cursor position.
//
// The vtkBrokenLineWidget has several methods that can be used in conjunction with
// other VTK objects. The Set/GetResolution() methods control the number of
// subdivisions of the broken line; the GetPolyData() method can be used to get the
// polygonal representation and can be used for things like seeding
// streamlines or probing other data sets. Typical usage of the widget is to
// make use of the StartInteractionEvent, InteractionEvent, and
// EndInteractionEvent events. The InteractionEvent is called on mouse motion;
// the other two events are called on button down and button up (either left or
// right button).
//
// Some additional features of this class include the ability to control the
// properties of the widget. You can set the properties of the selected and
// unselected representations of the broken line. For example, you can set the
// property for the handles and broken line. In addition there are methods to
// constrain the broken line so that it is aligned with a plane.  Note that a simple
// ruler widget can be derived by setting the resolution to 1, the number of
// handles to 2, and calling the GetSummedLength method!

// .SECTION Thanks
// Thanks to Philippe Pebay for developing and contributing this class.

// .SECTION Caveats
// Note that handles and line can be picked even when they are "behind" other
// actors.  This is an intended feature and not a bug.

// .SECTION See Also
// vtk3DWidget vtkBoxWidget vtkLineWidget vtkPointWidget vtkSphereWidget
// vtkImagePlaneWidget vtkImplicitPlaneWidget vtkPlaneWidget


#ifndef __vtkBrokenLineWidget_h
#define __vtkBrokenLineWidget_h

#include "vtk3DWidget.h"

class vtkActor;
class vtkCellPicker;
class vtkPlaneSource;
class vtkPoints;
class vtkPolyData;
class vtkProp;
class vtkProperty;
class vtkSphereSource;
class vtkTransform;

#define VTK_PROJECTION_YZ 0
#define VTK_PROJECTION_XZ 1
#define VTK_PROJECTION_XY 2
#define VTK_PROJECTION_OBLIQUE 3

class VTK_WIDGETS_EXPORT vtkBrokenLineWidget : public vtk3DWidget
{
public:
  // Description:
  // Instantiate the object.
  static vtkBrokenLineWidget *New();

  vtkTypeMacro(vtkBrokenLineWidget,vtk3DWidget);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Methods that satisfy the superclass' API.
  virtual void SetEnabled(int);
  virtual void PlaceWidget(double bounds[6]);
  void PlaceWidget()
    {this->Superclass::PlaceWidget();}
  void PlaceWidget(double xmin, double xmax, double ymin, double ymax, 
                   double zmin, double zmax)
    {this->Superclass::PlaceWidget(xmin,xmax,ymin,ymax,zmin,zmax);}

  // Description:
  // Force the broken line widget to be projected onto one of the orthogonal planes.
  // Remember that when the state changes, a ModifiedEvent is invoked.
  // This can be used to snap the broken line to the plane if it is orginally
  // not aligned.  The normal in SetProjectionNormal is 0,1,2 for YZ,XZ,XY
  // planes respectively and 3 for arbitrary oblique planes when the widget
  // is tied to a vtkPlaneSource.
  vtkSetMacro(ProjectToPlane,int);
  vtkGetMacro(ProjectToPlane,int);
  vtkBooleanMacro(ProjectToPlane,int);

  // Description:
  // Set up a reference to a vtkPlaneSource that could be from another widget
  // object, e.g. a vtkPolyDataSourceWidget.
  void SetPlaneSource(vtkPlaneSource* plane);

  vtkSetClampMacro(ProjectionNormal,int,VTK_PROJECTION_YZ,VTK_PROJECTION_OBLIQUE);
  vtkGetMacro(ProjectionNormal,int);
  void SetProjectionNormalToXAxes()
    { this->SetProjectionNormal(0); }
  void SetProjectionNormalToYAxes()
    { this->SetProjectionNormal(1); }
  void SetProjectionNormalToZAxes()
    { this->SetProjectionNormal(2); }
  void SetProjectionNormalToOblique()
    { this->SetProjectionNormal(3); }

  // Description:
  // Set the position of broken line handles and points in terms of a plane's
  // position. i.e., if ProjectionNormal is 0, all of the x-coordinate
  // values of the points are set to position. Any value can be passed (and is
  // ignored) to update the broken line points when Projection normal is set to 3
  // for arbritrary plane orientations.
  void SetProjectionPosition(double position);
  vtkGetMacro(ProjectionPosition, double);

  // Description:
  // Grab the polydata (including points) that defines the broken line.  The
  // polydata consists of points and line segments numbering Resolution + 1
  // and Resoltuion, respectively. Points are guaranteed to be up-to-date when
  // either the InteractionEvent or  EndInteraction events are invoked. The
  // user provides the vtkPolyData and the points and polyline are added to it.
  void GetPolyData(vtkPolyData *pd);

  // Description:
  // Set/Get the handle properties (the spheres are the handles). The
  // properties of the handles when selected and unselected can be manipulated.
  virtual void SetHandleProperty(vtkProperty*);
  vtkGetObjectMacro(HandleProperty, vtkProperty);
  virtual void SetSelectedHandleProperty(vtkProperty*);
  vtkGetObjectMacro(SelectedHandleProperty, vtkProperty);

  // Description:
  // Set/Get the line properties. The properties of the line when selected
  // and unselected can be manipulated.
  virtual void SetLineProperty(vtkProperty*);
  vtkGetObjectMacro(LineProperty, vtkProperty);
  virtual void SetSelectedLineProperty(vtkProperty*);
  vtkGetObjectMacro(SelectedLineProperty, vtkProperty);

  // Description:
  // Set/Get the number of handles for this widget.
  virtual void SetNumberOfHandles(int npts);
  vtkGetMacro(NumberOfHandles, int);

  // Description:
  // Set/Get the number of line segments representing the broken line for
  // this widget.
  // A broken line with resolution greater than the number of handles is useful when
  // points along the line are desired; e.g., generating a rake of streamlines.
  void SetResolution(int resolution);
  vtkGetMacro(Resolution,int);

  // Description:
  // Set/Get the position of the broken line handles. Call GetNumberOfHandles
  // to determine the valid range of handle indices.
  void SetHandlePosition(int handle, double x, double y, double z);
  void SetHandlePosition(int handle, double xyz[3]);
  void GetHandlePosition(int handle, double xyz[3]);
  double* GetHandlePosition(int handle);

  // Description:
  // Control whether the broken line is open or closed. A closed broken line forms
  // a continuous loop: the first and last points are the same, and
  // derivatives are continuous.  A minimum of 3 handles are required to
  // form a closed loop.  This method enforces consistency with
  // user supplied subclasses of vtkBrokenLine.
  void SetClosed(int closed);
  vtkGetMacro(Closed,int);
  vtkBooleanMacro(Closed,int);

  // Description:
  // Convenience method to determine whether the broken line is
  // closed in a geometric sense.  The widget may be set "closed" but still
  // be geometrically open (e.g., a straight line).
  int IsClosed();

  // Description:
  // Get the summed lengths of the individual straight line segments.
  double GetSummedLength();

  // Description:
  // Convenience method to allocate and set the handles from a vtkPoints
  // instance.  If the first and last points are the same, the broken line sets
  // Closed to the on state and disregards the last point, otherwise Closed
  // remains unchanged.
  void InitializeHandles(vtkPoints* points);

  // Description:
  // Turn on / off event processing for this widget. If off, the widget will
  // not respond to user interaction
  vtkSetClampMacro(ProcessEvents, int, 0, 1);
  vtkGetMacro(ProcessEvents, int);
  vtkBooleanMacro( ProcessEvents, int );

protected:
  vtkBrokenLineWidget();
  ~vtkBrokenLineWidget();

//BTX - manage the state of the widget
  int State;
  enum WidgetState
  {
    Start=0,
    Moving,
    Scaling,
    Spinning,
    Inserting,
    Erasing,
    Outside
  };
//ETX

  //handles the events
  static void ProcessEventsHandler(vtkObject* object,
                                   unsigned long event,
                                   void* clientdata,
                                   void* calldata);

  // ProcessEventsHandler() dispatches to these methods.
  void OnLeftButtonDown();
  void OnLeftButtonUp();
  void OnMiddleButtonDown();
  void OnMiddleButtonUp();
  void OnRightButtonDown();
  void OnRightButtonUp();
  void OnMouseMove();

  // Controlling vars
  int   ProjectionNormal;
  double ProjectionPosition;
  int   ProjectToPlane;
  vtkPlaneSource* PlaneSource;

  // Projection capabilities
  void ProjectPointsToPlane();
  void ProjectPointsToOrthoPlane();
  void ProjectPointsToObliquePlane();

  // The broken line
  int NumberOfHandles;
  int Closed;
  void BuildRepresentation();
  
  // The line segments
  vtkActor           *LineActor;
  void HighlightLine(int highlight);
  int Resolution;

  // Glyphs representing hot spots (e.g., handles)
  vtkActor          **Handle;
  vtkSphereSource   **HandleGeometry;
  void Initialize();
  int  HighlightHandle(vtkProp *prop); //returns handle index or -1 on fail
  virtual void SizeHandles();
  void InsertHandleOnLine(double* pos);
  void EraseHandle(const int&);

  // Do the picking
  vtkCellPicker *HandlePicker;
  vtkCellPicker *LinePicker;
  vtkActor *CurrentHandle;
  int CurrentHandleIndex;

  // Methods to manipulate the broken line.
  void MovePoint(double *p1, double *p2);
  void Scale(double *p1, double *p2, int X, int Y);
  void Translate(double *p1, double *p2);
  void Spin(double *p1, double *p2, double *vpn);

  // Transform the control points (used for spinning)
  vtkTransform *Transform;

  // Properties used to control the appearance of selected objects and
  // the manipulator in general.
  vtkProperty *HandleProperty;
  vtkProperty *SelectedHandleProperty;
  vtkProperty *LineProperty;
  vtkProperty *SelectedLineProperty;
  void CreateDefaultProperties();

  // For efficient spinning
  double Centroid[3];
  void CalculateCentroid();
  int  ProcessEvents;

private:
  vtkBrokenLineWidget(const vtkBrokenLineWidget&);  //Not implemented
  void operator=(const vtkBrokenLineWidget&);  //Not implemented
};

#endif
