/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkPlot3D.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

// .NAME vtkPlot3D - Abstract class for 3D plots.
//
// .SECTION Description
// The base class for all plot types used in vtkChart derived charts.
//
// .SECTION See Also
// vtkPlot3DPoints vtkPlot3DLine vtkPlot3DBar vtkChart vtkChartXY

#ifndef __vtkPlot3D_h
#define __vtkPlot3D_h

#include "vtkChartsCoreModule.h" // For export macro
#include "vtkContextItem.h"
#include "vtkNew.h"              // Needed to hold vtkNew ivars
#include "vtkSmartPointer.h"     // Needed to hold SP ivars
#include "vtkVector.h"           // For Points ivar
#include <vector>                // For ivars
using std::vector;

class vtkChartXYZ;
class vtkTable;
class vtkUnsignedCharArray;
class vtkPen;

class VTKCHARTSCORE_EXPORT vtkPlot3D : public vtkContextItem
{
public:
  vtkTypeMacro(vtkPlot3D, vtkContextItem);
  virtual void PrintSelf(ostream &os, vtkIndent indent);

  // Description:
  // Set/get the vtkPen object that controls how this plot draws (out)lines.
  void SetPen(vtkPen *pen);
  vtkPen* GetPen();

  // Description:
  // Set the input to the plot.
  virtual void SetInputData(vtkTable *input);
  virtual void SetInputData(vtkTable *input, const vtkStdString &xName,
                            const vtkStdString &yName,
                            const vtkStdString &zName);
  virtual void SetInputData(vtkTable *input, const vtkStdString &xName,
                            const vtkStdString &yName,
                            const vtkStdString &zName,
                            const vtkStdString &colorName);
  virtual void SetInputData(vtkTable *input, vtkIdType xColumn,
                            vtkIdType yColumn, vtkIdType zColumn);

  // Description:
  // Get all the data points within this plot.
  vector<vtkVector3f> GetPoints();

  // Description:
  // Get the data points from this plot that currently fall within the axes
  // of the chart.
  vector<vtkVector3f> GetPointsThatSurviveClipping();

  // Description:
  // Get/set the chart for this plot.
  vtkGetObjectMacro(Chart, vtkChartXYZ);
  virtual void SetChart(vtkChartXYZ* chart);

  // Description:
  // Get the label for the X axis.
  std::string GetXAxisLabel();

  // Description:
  // Get the label for the Y axis.
  std::string GetYAxisLabel();

  // Description:
  // Get the label for the Z axis.
  std::string GetZAxisLabel();

  // Description:
  // Get the bounding cube surrounding the currently rendered data points.
  vector<vtkVector3f> GetDataBounds() { return this->DataBounds; }

//BTX
protected:
  vtkPlot3D();
  ~vtkPlot3D();

  // Description:
  // Determine what data points fall within the bounds of the chart axes.
  virtual void UpdateClippedPoints();

  // Description:
  // Generate a bounding cube for our data.
  virtual void ComputeDataBounds();

  // Description:
  // This object stores the vtkPen that controls how the plot is drawn.
  vtkSmartPointer<vtkPen> Pen;

  // Description:
  // This array assigns a color to each datum in the plot.
  vtkNew<vtkUnsignedCharArray> Colors;

  // Description:
  // This array assigns a color to each data point which is currently
  // rendered within the axes.
  vtkNew<vtkUnsignedCharArray> ClippedColors;

  // Description:
  // Number of components in our color vectors.  This value is initialized
  // to zero.  It's typically set to 3 or 4 if the points are to be colored.
  int NumberOfComponents;

  // Description:
  // The label for the X Axis.
  std::string XAxisLabel;

  // Description:
  // The label for the Y Axis.
  std::string YAxisLabel;

  // Description:
  // The label for the Z Axis.
  std::string ZAxisLabel;

  // Description:
  // The data points read in during SetInputData().
  vector<vtkVector3f> Points;

  // Description:
  // The subset of our data points that fall within the axes.  These are
  // the only data points that are rendered.
  vector<vtkVector3f> PointsThatSurviveClipping;

  // Description:
  // When the points were last built.
  vtkTimeStamp PointsBuildTime;

  // Description:
  // The chart containing this plot.
  vtkChartXYZ* Chart;

  // Description:
  // A bounding cube surrounding the currently rendered data points.
  vector<vtkVector3f> DataBounds;

private:
  vtkPlot3D(const vtkPlot3D &); // Not implemented.
  void operator=(const vtkPlot3D &); // Not implemented.

//ETX
};

#endif //__vtkPlot3D_h
