/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkInteractiveChartXYZ.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

// .NAME vtkInteractiveChartXYZ - Factory class for drawing 3D XYZ charts.
//
// .SECTION Description

#ifndef __vtkInteractiveChartXYZ_h
#define __vtkInteractiveChartXYZ_h

#include "vtkChartsCoreModule.h" // For export macro
#include "vtkChartXYZ.h"
#include "vtkRect.h"         // For vtkRectf ivars
#include "vtkNew.h"          // For ivars
#include "vtkSmartPointer.h" // For ivars

class vtkAnnotationLink;
class vtkAxis;
class vtkContextMouseEvent;
class vtkPlot;
class vtkTable;
class vtkTransform;
class vtkPen;

class VTKCHARTSCORE_EXPORT vtkInteractiveChartXYZ : public vtkChartXYZ
{
public:
  vtkTypeMacro(vtkInteractiveChartXYZ, vtkChartXYZ);
  virtual void PrintSelf(ostream &os, vtkIndent indent);

  static vtkInteractiveChartXYZ * New();

  // Description:
  // Update any data as necessary before drawing the chart.
  void Update();

  // Description:
  // Paint event for the chart, called whenever the chart needs to be drawn
  virtual bool Paint(vtkContext2D *painter);

  // Description:
  // Set the input for the chart, this should be done in the plot, but keeping
  // things simple while I get everything working...
  virtual void SetInput(vtkTable *input, const vtkStdString &x,
                        const vtkStdString &y, const vtkStdString &z);
  virtual void SetInput(vtkTable *input, const vtkStdString &x,
                        const vtkStdString &y, const vtkStdString &z,
                        const vtkStdString &r,  const vtkStdString &g,
                        const vtkStdString &b);
  virtual void SetInput(vtkTable *input, const vtkStdString &x,
                        const vtkStdString &y, const vtkStdString &z,
                        const vtkStdString &r,  const vtkStdString &g,
                        const vtkStdString &b, const vtkStdString &a);

  //BTX
  // Description:
  // Returns true if the transform is interactive, false otherwise.
  virtual bool Hit(const vtkContextMouseEvent &mouse);

  // Description:
  // Mouse press event. Keep track of zoom anchor position.
  virtual bool MouseButtonPressEvent(const vtkContextMouseEvent &mouse);

  // Description:
  // Mouse move event. Perform pan or zoom as specified by the mouse bindings.
  virtual bool MouseMoveEvent(const vtkContextMouseEvent &mouse);

  // Description:
  // Mouse wheel event. Perform pan or zoom as specified by mouse bindings.
  virtual bool MouseWheelEvent(const vtkContextMouseEvent &mouse, int delta);
  //ETX

protected:
  vtkInteractiveChartXYZ();
  ~vtkInteractiveChartXYZ();
  virtual void CalculateTransforms();
  bool Rotate(const vtkContextMouseEvent &mouse);
  bool Pan(const vtkContextMouseEvent &mouse);
  bool Zoom(const vtkContextMouseEvent &mouse);
  bool Spin(const vtkContextMouseEvent &mouse);

  vtkNew<vtkTransform> Translation;
  vtkNew<vtkTransform> Scale;
  unsigned char *Colors;
  int NumberOfComponents;

private:
  vtkInteractiveChartXYZ(const vtkInteractiveChartXYZ &);    // Not implemented.
  void operator=(const vtkInteractiveChartXYZ &); // Not implemented.
};

#endif
