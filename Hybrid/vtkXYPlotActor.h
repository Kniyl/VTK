/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkXYPlotActor.h
  Language:  C++
  Date:      $Date$
  Version:   $Revision$
  Thanks:    Thanks to Kitware & RPI/SCOREC who supported the development
             of this class.

Copyright (c) 1993-2001 Ken Martin, Will Schroeder, Bill Lorensen 
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.

 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

 * Neither name of Ken Martin, Will Schroeder, or Bill Lorensen nor the names
   of any contributors may be used to endorse or promote products derived
   from this software without specific prior written permission.

 * Modified source versions must be plainly marked as such, and must not be
   misrepresented as being the original software.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS ``AS IS''
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=========================================================================*/
// .NAME vtkXYPlotActor - generate an x-y plot from input dataset(s) or field data
// .SECTION Description
// vtkXYPlotActor creates an x-y plot of data from one or more input data
// sets or field data. The class plots dataset scalar values (y-axis) against
// the points (x-axis). The x-axis values are generated by taking the point
// ids, computing a cumulative arc length, or a normalized arc length. More
// than one input data set can be specified to generate multiple plots.
// Alternatively, if field data is supplied as input, the class plots one
// component against another. (The user must specify which component to use
// as the x-axis and which for the y-axis.)
//
// To use this class to plot dataset(s), you must specify one or more
// input datasets containing scalar and point data.  You'll probably also
// want to invoke a method to control how the point coordinates are converted
// into x values (by default point ids are used).
//
// To use this class to plot field data, you must specify one or more input
// data objects with its associated field data. You'll also want to specify
// which component to use as the x-axis and which to use as the y-axis.
// Note that when plotting field data, the x and y values are used directly
// (i.e., there are no options to normalize the components).
//
// Once you've set up the plot, you'll want to position it.  The
// PositionCoordinate defines the lower-left location of the x-y plot
// (specified in normalized viewport coordinates) and the Position2Coordinate
// define the upper-right corner. (Note: the Position2Coordinate is relative
// to PositionCoordinate, so you can move the vtkXYPlotActor around the
// viewport by setting just the PositionCoordinate.) The combination of the
// two position coordinates specifies a rectangle in which the plot will lie.
//
// Optional features include the ability to specify axes labels, label
// format, plot title, and control font properties and type. You can also
// manually specify the x and y plot ranges (by default they are computed
// automatically). The Border instance variable is used to create space 
// between the boundary of the plot window (specified by PositionCoordinate
// and Position2Coordinate) and the plot itself.
//
// There are several advanced features as well. You can assign per curve 
// properties (such as color and a plot symbol). (Note that each input 
// dataset and/or data object creates a single curve.) Another option is to
// add a plot legend that graphically indicates the correspondance between
// the curve, curve symbols, and the data source. You can also exchange the
// x and y axes if you prefer you plot orientation that way.

// .SECTION Caveats
// If you are interested in plotting something other than scalar data, you
// can use the vtk data shuffling filters (e.g., 
// vtkAttributeDataToFieldDataFilter snd vtkFieldDataToAttributeDataFilter) 
// to convert the data into scalar data and/or points.

// .SECTION See Also
// vtkActor2D vtkTextMapper vtkScalarBarActor vtkAxisActor2D vtkCubeAxesActor
// vtkAttributeDataToFieldDataFilter vtkFieldDataToAttributeDataFilter

#ifndef __vtkXYPlotActor_h
#define __vtkXYPlotActor_h

#define VTK_XYPLOT_INDEX                 0
#define VTK_XYPLOT_ARC_LENGTH            1
#define VTK_XYPLOT_NORMALIZED_ARC_LENGTH 2
#define VTK_XYPLOT_VALUE                 3

#define VTK_XYPLOT_ROW 0
#define VTK_XYPLOT_COLUMN 1

#include "vtkAxisActor2D.h"
class vtkDataSetCollection;
class vtkDataObjectCollection;
class vtkGlyphSource2D;
class vtkGlyph2D;
class vtkLegendBoxActor;
class vtkAppendPolyData;
class vtkPlanes;
class vtkIntArray;

class VTK_HYBRID_EXPORT vtkXYPlotActor : public vtkActor2D
{
public:
  vtkTypeRevisionMacro(vtkXYPlotActor,vtkActor2D);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Instantiate object with autorange computation; bold, italic, and shadows
  // on; arial font family; the number of labels set to 5 for the x and y
  // axes; a label format of "%-#6.3g"; and x coordinates computed from point
  // ids.
  static vtkXYPlotActor *New();

  //---Data Set Input----------------------------------------------------------
  // The following methods are used to plot input datasets. Datasets
  // will be plotted if set as input; otherwise the input data objects
  // will be plotted (if defined).
  
  // Description:
  // Add a dataset to the list of data to append. The array name specifies
  // which point array to plot.  If the array name is NULL, then the default
  // scalars are used.  The array can have multiple components, but only the
  // first component is ploted.
  void AddInput(vtkDataSet *in, const char* arrayName, int component);
  void AddInput(vtkDataSet *in) {this->AddInput(in, NULL, 0);}

  // Description:
  // Remove a dataset from the list of data to append.
  void RemoveInput(vtkDataSet *in, const char* arrayName, int component);
  void RemoveInput(vtkDataSet *in) {this->RemoveInput(in, NULL, 0);}

  // Description:
  // This removes all of the data set inputs, 
  // but does not change the data object inputs.
  void RemoveAllInputs();

  // Description:
  // Return the list of inputs to this filter.
  vtkDataSetCollection *GetInputList() {return this->InputList;}

  // Description:
  // If plotting points by value, which component to use to determine the
  // value. This sets a value per each input dataset (i.e., the ith dataset).
  void SetPointComponent(int i, int comp);
  int GetPointComponent(int i);
  //---end Data Set Input-----------------------------------------------------

  // Description:
  // Specify how the independent (x) variable is computed from the points.
  // The independent variable can be the scalar/point index (i.e., point id),
  // the accumulated arc length along the points, the normalized arc length,
  // or by component value. If plotting datasets (e.g., points), the value
  // that is used is specified by the PointComponent ivar.  (Note: these
  // methods also control how field data is plotted. Field data is usually
  // plotted by value or index, if plotting length 1-dimensional length
  // measures are used.)
  vtkSetClampMacro(XValues,int,VTK_XYPLOT_INDEX,VTK_XYPLOT_VALUE);
  vtkGetMacro(XValues,int);
  void SetXValuesToIndex(){this->SetXValues(VTK_XYPLOT_INDEX);};
  void SetXValuesToArcLength() {this->SetXValues(VTK_XYPLOT_ARC_LENGTH);};
  void SetXValuesToNormalizedArcLength()
    {this->SetXValues(VTK_XYPLOT_NORMALIZED_ARC_LENGTH);};
  void SetXValuesToValue() {this->SetXValues(VTK_XYPLOT_VALUE);};
  const char *GetXValuesAsString();

  //---Data Object Input------------------------------------------------------
  // The following methods are used to plot input data objects. Datasets will
  // be plotted in preference to data objects if set as input; otherwise the
  // input data objects will be plotted (if defined).
  
  // Description:
  // Add a dataset to the list of data to append.
  void AddDataObjectInput(vtkDataObject *in);

  // Description:
  // Remove a dataset from the list of data to append.
  void RemoveDataObjectInput(vtkDataObject *in);

  // Description:
  // Return the list of inputs to this filter.
  vtkDataObjectCollection *GetDataObjectInputList() 
    {return this->DataObjectInputList;}

  // Description:
  // Indicate whether to plot rows or columns. If plotting rows, then
  // the dependent variables is taken from a specified row,
  // versus rows (y). 
  vtkSetClampMacro(DataObjectPlotMode,int,VTK_XYPLOT_COLUMN,VTK_XYPLOT_ROW);
  vtkGetMacro(DataObjectPlotMode,int);
  void SetDataObjectPlotModeToRows()
    {this->SetDataObjectPlotMode(VTK_XYPLOT_ROW);}
  void SetDataObjectPlotModeToColumns()
    {this->SetDataObjectPlotMode(VTK_XYPLOT_COLUMN);}
  const char *GetDataObjectPlotModeAsString();

  // Description:
  // Specify which component of the input data object to use as the
  // independent variable for the ith input data object. (This ivar is
  // ignored if plotting the index.) Note that the value is interpreted
  // differently depending on DataObjectPlotMode. If the mode is Rows, then
  // the value of DataObjectXComponent is the row number; otherwise it's the
  // column number.
  void SetDataObjectXComponent(int i, int comp);
  int GetDataObjectXComponent(int i);

  // Description:
  // Specify which component of the input data object to use as the
  // dependent variable for the ith input data object. (This ivar is
  // ignored if plotting the index.) Note that the value is interpreted
  // differently depending on DataObjectPlotMode. If the mode is Rows, then
  // the value of DataObjectYComponent is the row number; otherwise it's the
  // column number.
  void SetDataObjectYComponent(int i, int comp);
  int GetDataObjectYComponent(int i);
  //---end Data Object Input--------------------------------------------------

  //---Per Curve Properties---------------------------------------------------
  // The following methods are used to set properties on each curve that is
  // plotted. Each input dataset (or data object) results in one curve. The
  // methods that follow have an index i that corresponds to the input dataset
  // or data object. 
  void SetPlotColor(int i, float r, float g, float b);
  void SetPlotColor(int i, const float color[3]) {
    this->SetPlotColor(i, color[0], color[1], color[2]); };
  float *GetPlotColor(int i);
  void SetPlotSymbol(int i,vtkPolyData *input);
  vtkPolyData *GetPlotSymbol(int i);
  void SetPlotLabel(int i, const char *label);
  const char *GetPlotLabel(int i);

  // Allow per-curve specification of line and point rendering.  These override
  // global settings PlotPoints and PlotLines.  If not on, the default behavior
  // is governed by PlotPoints and PlotLines ivars.
  vtkGetMacro(PlotCurvePoints, int);
  vtkSetMacro(PlotCurvePoints, int);
  vtkBooleanMacro(PlotCurvePoints, int);

  vtkGetMacro(PlotCurveLines, int);
  vtkSetMacro(PlotCurveLines, int);
  vtkBooleanMacro(PlotCurveLines, int);

  void SetPlotLines(int i, int);
  int GetPlotLines(int i);

  void SetPlotPoints(int i, int);
  int GetPlotPoints(int i);
  //---end Per Curve Properties-----------------------------------------------

  // Description:
  // Enable/Disable exchange of the x-y axes (i.e., what was x becomes y, and
  // vice-versa). Exchanging axes affects the labeling as well.
  vtkSetMacro(ExchangeAxes, int);
  vtkGetMacro(ExchangeAxes, int);
  vtkBooleanMacro(ExchangeAxes, int);

  // Description:
  // Normally the x-axis is plotted from minimum to maximum. Setting this instance
  // variable causes the x-axis to be plotted from maximum to minimum. Note that
  // boolean always applies to the x-axis even if ExchangeAxes is set.
  vtkSetMacro(ReverseXAxis, int);
  vtkGetMacro(ReverseXAxis, int);
  vtkBooleanMacro(ReverseXAxis, int);

  // Description:
  // Normally the y-axis is plotted from minimum to maximum. Setting this instance
  // variable causes the y-axis to be plotted from maximum to minimum. Note that
  // boolean always applies to the y-axis even if ExchangeAxes is set.
  vtkSetMacro(ReverseYAxis, int);
  vtkGetMacro(ReverseYAxis, int);
  vtkBooleanMacro(ReverseYAxis, int);

  // Description:
  // Retrieve handles to the legend box and glyph source. This is useful
  // if you would like to change the default behavior of the legend box
  // or glyph source. For example, the default glyph can be changed from
  // a line to a vertex plus line, etc.)
  vtkLegendBoxActor *GetLegendBoxActor()
    {return this->LegendActor;}
  vtkGlyphSource2D *GetGlyphSource()
    {return this->GlyphSource;}

  // Description:
  // Set/Get the title of the x-y plot, and the title along the 
  // x and y axes.
  vtkSetStringMacro(Title);
  vtkGetStringMacro(Title);
  vtkSetStringMacro(XTitle);
  vtkGetStringMacro(XTitle);
  vtkSetStringMacro(YTitle);
  vtkGetStringMacro(YTitle);

  // Description:
  // Set the plot range (range of independent and dependent variables)
  // to plot. Data outside of the range will be clipped. If the plot
  // range of either the x or y variables is set to (v1,v2), where
  // v1 == v2, then the range will be computed automatically. Note that
  // the x-range values should be consistent with the way the independent
  // variable is created (via INDEX, DISTANCE, or ARC_LENGTH).
  vtkSetVector2Macro(XRange,float);
  vtkGetVectorMacro(XRange,float,2);
  vtkSetVector2Macro(YRange,float);
  vtkGetVectorMacro(YRange,float,2);
  void SetPlotRange(float xmin, float ymin, float xmax, float ymax)
    {this->SetXRange(xmin,xmax); this->SetYRange(ymin,ymax);}
  
  // Description:
  // Set/Get the number of annotation labels to show along the x and y axes.
  // This values is a suggestion: the number of labels may vary depending
  // on the particulars of the data. The convenience method 
  // SetNumberOfLables() sets the number of x and y labels to the same value.
  vtkSetClampMacro(NumberOfXLabels, int, 0, 50);
  vtkGetMacro(NumberOfXLabels, int);
  vtkSetClampMacro(NumberOfYLabels, int, 0, 50);
  vtkGetMacro(NumberOfYLabels, int);
  void SetNumberOfLabels(int num)
    {this->SetNumberOfXLabels(num); this->SetNumberOfYLabels(num);}
  
  // Description:
  // Enable/Disable the creation of a legend. If on, the legend labels will
  // be created automatically unless the per plot legend symbol has been
  // set.
  vtkSetMacro(Legend, int);
  vtkGetMacro(Legend, int);
  vtkBooleanMacro(Legend, int);

  // Description: 
  // Use these methods to control the position of the legend. The variables
  // LegendPosition and LegendPosition2 define the lower-left and upper-right
  // position of the legend. The coordinates are expressed as normalized
  // values with respect to the rectangle defined by PositionCoordinate and
  // Position2Coordinate. Note that LegendPosition2 is relative to
  // LegendPosition.
  vtkSetVector2Macro(LegendPosition,float);
  vtkGetVector2Macro(LegendPosition,float);
  vtkSetVector2Macro(LegendPosition2,float);
  vtkGetVector2Macro(LegendPosition2,float);
  
  // Description:
  // Enable/Disable bolding annotation text.
  vtkSetMacro(Bold, int);
  vtkGetMacro(Bold, int);
  vtkBooleanMacro(Bold, int);

  // Description:
  // Enable/Disable italicizing annotation text.
  vtkSetMacro(Italic, int);
  vtkGetMacro(Italic, int);
  vtkBooleanMacro(Italic, int);

  // Description:
  // Enable/Disable creating shadows on the annotation text. Shadows make 
  // the text easier to read.
  vtkSetMacro(Shadow, int);
  vtkGetMacro(Shadow, int);
  vtkBooleanMacro(Shadow, int);

  // Description:
  // Enable/Disable plotting of Log of x-values.
  vtkSetMacro(Logx, int);
  vtkGetMacro(Logx, int);
  vtkBooleanMacro(Logx, int);

  // Description:
  // Set/Get the font family for the annotation text. Three font types 
  // are available: Arial (VTK_ARIAL), Courier (VTK_COURIER), and 
  // Times (VTK_TIMES).
  vtkSetMacro(FontFamily, int);
  vtkGetMacro(FontFamily, int);
  void SetFontFamilyToArial() {this->SetFontFamily(VTK_ARIAL);};
  void SetFontFamilyToCourier() {this->SetFontFamily(VTK_COURIER);};
  void SetFontFamilyToTimes() {this->SetFontFamily(VTK_TIMES);};

  // Description:
  // Set/Get the format with which to print the labels on the scalar
  // bar.
  vtkSetStringMacro(LabelFormat);
  vtkGetStringMacro(LabelFormat);

  // Description:
  // Set/Get the spacing between the plot window and the plot. The value
  // is specified in pixels.
  vtkSetClampMacro(Border, int, 0, 50);
  vtkGetMacro(Border, int);

  // Description:
  // Set/Get whether the points are rendered.  The point size can be set in
  // the property object. This is a global flag which affects the plot only 
  // if per curve symbols are not defined.
  vtkGetMacro(PlotPoints, int);
  vtkSetMacro(PlotPoints, int);
  vtkBooleanMacro(PlotPoints, int);

  // Description:
  // Set/Get whether the lines are rendered.  The line width can be set in
  // the property object. 
  vtkGetMacro(PlotLines, int);
  vtkSetMacro(PlotLines, int);
  vtkBooleanMacro(PlotLines, int);
  
  // Description:
  // Set/Get the factor that controls how big glyphs are in the plot.
  // The number is expressed as a fraction of the length of the diagonal
  // of the plot bounding box.
  vtkSetClampMacro(GlyphSize, float, 0.0, 0.2);
  vtkGetMacro(GlyphSize, float);

  // Description:
  // Given a position within the viewport used by the plot, return the
  // the plot coordinates (XAxis value, YAxis value)
  void ViewportToPlotCoordinate(vtkViewport *viewport, float &u, float &v);

  // Description:
  // An alternate form of ViewportToPlotCoordinate() above. This method
  // inputs the viewport coordinate pair (defined by the ivar 
  // ViewportCoordinate)and then stores them in the ivar PlotCoordinate. 
  void ViewportToPlotCoordinate(vtkViewport *viewport);
  vtkSetVector2Macro(PlotCoordinate,float);
  vtkGetVector2Macro(PlotCoordinate,float);

  // Description:
  // Given a plot coordinate, return the viewpoint position
  void PlotToViewportCoordinate(vtkViewport *viewport, float &u, float &v);

  // Description:
  // An alternate form of PlotToViewportCoordinate() above. This method
  // inputs the plot coordinate pair (defined in the ivar PlotCoordinate)
  // and then stores them in the ivar ViewportCoordinate. (This method 
  // can be wrapped.)
  void PlotToViewportCoordinate(vtkViewport *viewport);
  vtkSetVector2Macro(ViewportCoordinate,float);
  vtkGetVector2Macro(ViewportCoordinate,float);

  // Description:
  // Is the specified viewport position within the plot area (as opposed to the
  // region used by the plot plus the labels)?
  int IsInPlot(vtkViewport *viewport, float u, float v);
  
  // Description:
  // Take into account the modified time of internal helper classes.
  unsigned long GetMTime();
  
//BTX  
  // Description:
  // WARNING: INTERNAL METHOD - NOT INTENDED FOR GENERAL USE
  // DO NOT USE THIS METHOD OUTSIDE OF THE RENDERING PROCESS.
  // Draw the x-y plot.
  int RenderOpaqueGeometry(vtkViewport*);
  int RenderOverlay(vtkViewport*);
  int RenderTranslucentGeometry(vtkViewport *) {return 0;}

  // Description:
  // Release any graphics resources that are being consumed by this actor.
  // The parameter window could be used to determine which graphic
  // resources to release.
  void ReleaseGraphicsResources(vtkWindow *);
//ETX  

protected:
  vtkXYPlotActor();
  ~vtkXYPlotActor();

  vtkDataSetCollection *InputList; //list of data sets to plot
  char** SelectedInputScalars; // list of data set arrays to plot
  vtkIntArray* SelectedInputScalarsComponent; // list of componenents
  vtkDataObjectCollection *DataObjectInputList; //list of data objects to plot
  char  *Title;
  char  *XTitle;
  char  *YTitle;
  int   XValues;
  int   NumberOfXLabels;
  int   NumberOfYLabels;
  int   Bold;
  int   Italic;
  int   Shadow;
  int   FontFamily;
  int   Logx;
  char  *LabelFormat;
  float XRange[2];
  float YRange[2];
  float XComputedRange[2];  //range actually used by plot
  float YComputedRange[2];  //range actually used by plot
  int Border;
  int PlotLines;
  int PlotPoints;
  int PlotCurveLines;
  int PlotCurvePoints;
  int ExchangeAxes;
  int ReverseXAxis;
  int ReverseYAxis;
  
  vtkTextMapper *TitleMapper;
  vtkActor2D    *TitleActor;

  vtkAxisActor2D      *XAxis;
  vtkAxisActor2D      *YAxis;

  float ViewportCoordinate[2];
  float PlotCoordinate[2];
  
  //Handle data objects and datasets
  int DataObjectPlotMode;
  vtkIntArray *XComponent;
  vtkIntArray *YComponent;
  vtkIntArray *LinesOn;
  vtkIntArray *PointsOn;

  //The data drawn within the axes. Each curve is one polydata.
  //color is controlled by scalar data. The curves are appended
  //together, possibly glyphed with point symbols.
  int NumberOfInputs;
  vtkPolyData             **PlotData; 
  vtkGlyph2D              **PlotGlyph;
  vtkAppendPolyData       **PlotAppend;
  vtkPolyDataMapper2D     **PlotMapper;
  vtkActor2D              **PlotActor;
  void                    InitializeEntries();
  
  // Legends and plot symbols. The legend also keeps track of
  // the symbols and such.
  int Legend;
  float LegendPosition[2];
  float LegendPosition2[2];
  vtkLegendBoxActor *LegendActor;
  vtkGlyphSource2D *GlyphSource;
  vtkPlanes *ClipPlanes;
  float GlyphSize;

  // Keep track of changes.
  int CachedSize[2];
  vtkTimeStamp  BuildTime;

  void ComputeXRange(float range[2], float *lengths);
  void ComputeYRange(float range[2]);
  void ComputeDORange(float xrange[2], float yrange[2], float *lengths);

  virtual void CreatePlotData(int *pos, int *pos2, float xRange[2], 
                              float yRange[2], float *norms, 
                              int numDS, int numDO);
  void PlaceAxes(vtkViewport *viewport, int *size, int pos[2], int pos2[2]);
  void GenerateClipPlanes(int *pos, int *pos2);
  float ComputeGlyphScale(int i, int *pos, int *pos2);
  void ClipPlotData(int *pos, int *pos2, vtkPolyData *pd);
  float *TransformPoint(int pos[2], int pos2[2], float x[3], float xNew[3]);
  
private:
  vtkXYPlotActor(const vtkXYPlotActor&);  // Not implemented.
  void operator=(const vtkXYPlotActor&);  // Not implemented.
};


#endif

