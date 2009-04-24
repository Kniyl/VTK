/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkRenderedGraphRepresentation.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
/*-------------------------------------------------------------------------
  Copyright 2008 Sandia Corporation.
  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
  the U.S. Government retains certain rights in this software.
-------------------------------------------------------------------------*/
// .NAME vtkRenderedGraphRepresentation - 
//
// .SECTION Description

#ifndef __vtkRenderedGraphRepresentation_h
#define __vtkRenderedGraphRepresentation_h

#include "vtkRenderedRepresentation.h"
#include "vtkSmartPointer.h" // for SP ivars

class vtkActor;
class vtkApplyColors;
class vtkArrayMap;
class vtkEdgeCenters;
class vtkEdgeLayout;
class vtkEdgeLayoutStrategy;
class vtkGraphLayout;
class vtkGraphLayoutStrategy;
class vtkGraphToGlyphs;
class vtkGraphToPoints;
class vtkGraphToPolyData;
class vtkInformation;
class vtkInformationVector;
class vtkLookupTable;
class vtkPerturbCoincidentVertices;
class vtkPolyData;
class vtkPolyDataMapper;
class vtkRenderView;
class vtkScalarBarWidget;
class vtkScalarsToColors;
class vtkTextProperty;
class vtkVertexDegree;
class vtkView;
class vtkViewTheme;

class VTK_VIEWS_EXPORT vtkRenderedGraphRepresentation : public vtkRenderedRepresentation
{
public:
  static vtkRenderedGraphRepresentation* New();
  vtkTypeRevisionMacro(vtkRenderedGraphRepresentation, vtkRenderedRepresentation);
  void PrintSelf(ostream& os, vtkIndent indent);

  virtual void SetupRenderWindow(vtkRenderWindow* win);

  virtual void SetVertexLabelArrayName(const char* name);
  virtual const char* GetVertexLabelArrayName();
  virtual void SetVertexLabelPriorityArrayName(const char* name);
  virtual const char* GetVertexLabelPriorityArrayName();
  virtual void SetVertexLabelVisibility(bool b);
  virtual bool GetVertexLabelVisibility();
  vtkBooleanMacro(VertexLabelVisibility, bool);
  virtual void SetVertexLabelTextProperty(vtkTextProperty* p);
  virtual vtkTextProperty* GetVertexLabelTextProperty();

  virtual void SetEdgeLabelArrayName(const char* name);
  virtual const char* GetEdgeLabelArrayName();
  virtual void SetEdgeLabelPriorityArrayName(const char* name);
  virtual const char* GetEdgeLabelPriorityArrayName();
  virtual void SetEdgeLabelVisibility(bool b);
  virtual bool GetEdgeLabelVisibility();
  vtkBooleanMacro(EdgeLabelVisibility, bool);
  virtual void SetEdgeLabelTextProperty(vtkTextProperty* p);
  virtual vtkTextProperty* GetEdgeLabelTextProperty();

  virtual void SetVertexIconArrayName(const char* name);
  virtual const char* GetVertexIconArrayName();
  virtual void SetVertexIconPriorityArrayName(const char* name);
  virtual const char* GetVertexIconPriorityArrayName();
  virtual void SetVertexIconVisibility(bool b);
  virtual bool GetVertexIconVisibility();
  vtkBooleanMacro(VertexIconVisibility, bool);
  virtual void AddVertexIconType(const char* name, int type);
  virtual void ClearVertexIconTypes();
  virtual void SetUseVertexIconTypeMap(bool b);
  virtual bool GetUseVertexIconTypeMap();
  vtkBooleanMacro(UseVertexIconTypeMap, bool);
  virtual void SetVertexIconAlignment(int align);
  virtual int GetVertexIconAlignment();

  virtual void SetEdgeIconArrayName(const char* name);
  virtual const char* GetEdgeIconArrayName();
  virtual void SetEdgeIconPriorityArrayName(const char* name);
  virtual const char* GetEdgeIconPriorityArrayName();
  virtual void SetEdgeIconVisibility(bool b);
  virtual bool GetEdgeIconVisibility();
  vtkBooleanMacro(EdgeIconVisibility, bool);
  virtual void AddEdgeIconType(const char* name, int type);
  virtual void ClearEdgeIconTypes();
  virtual void SetUseEdgeIconTypeMap(bool b);
  virtual bool GetUseEdgeIconTypeMap();
  vtkBooleanMacro(UseEdgeIconTypeMap, bool);
  virtual void SetEdgeIconAlignment(int align);
  virtual int GetEdgeIconAlignment();

  virtual void SetColorVerticesByArray(bool b);
  virtual bool GetColorVerticesByArray();
  vtkBooleanMacro(ColorVerticesByArray, bool);
  virtual void SetVertexColorArrayName(const char* name);
  virtual const char* GetVertexColorArrayName();

  virtual void SetColorEdgesByArray(bool b);
  virtual bool GetColorEdgesByArray();
  vtkBooleanMacro(ColorEdgesByArray, bool);
  virtual void SetEdgeColorArrayName(const char* name);
  virtual const char* GetEdgeColorArrayName();

  virtual void SetEnableVerticesByArray(bool b);
  virtual bool GetEnableVerticesByArray();
  vtkBooleanMacro(EnableVerticesByArray, bool);
  virtual void SetEnabledVerticesArrayName(const char* name);
  virtual const char* GetEnabledVerticesArrayName();

  virtual void SetEnableEdgesByArray(bool b);
  virtual bool GetEnableEdgesByArray();
  vtkBooleanMacro(EnableEdgesByArray, bool);
  virtual void SetEnabledEdgesArrayName(const char* name);
  virtual const char* GetEnabledEdgesArrayName();

  virtual void SetEdgeVisibility(bool b);
  virtual bool GetEdgeVisibility();
  vtkBooleanMacro(EdgeVisibility, bool);

  // Description:
  // Set/get the graph layout strategy.
  virtual void SetLayoutStrategy(vtkGraphLayoutStrategy* strategy);
  virtual vtkGraphLayoutStrategy* GetLayoutStrategy();

  // Description:
  // Get/set the layout strategy by name.
  virtual void SetLayoutStrategy(const char* name);
  vtkGetStringMacro(LayoutStrategyName);

  // Description:
  // Set predefined layout strategies.
  void SetLayoutStrategyToRandom()
    { this->SetLayoutStrategy("Random"); }
  void SetLayoutStrategyToForceDirected()
    { this->SetLayoutStrategy("Force Directed"); }
  void SetLayoutStrategyToSimple2D()
    { this->SetLayoutStrategy("Simple 2D"); }
  void SetLayoutStrategyToClustering2D()
    { this->SetLayoutStrategy("Clustering 2D"); }
  void SetLayoutStrategyToCommunity2D()
    { this->SetLayoutStrategy("Community 2D"); }
  void SetLayoutStrategyToFast2D()
    { this->SetLayoutStrategy("Fast 2D"); }
  void SetLayoutStrategyToPassThrough()
    { this->SetLayoutStrategy("Pass Through"); }
  void SetLayoutStrategyToCircular()
    { this->SetLayoutStrategy("Circular"); }

  // Description:
  // Set the layout strategy to use coordinates from arrays.
  // The x array must be specified. The y and z arrays are optional.
  virtual void SetLayoutStrategyToAssignCoordinates(
    const char* xarr, const char* yarr = 0, const char* zarr = 0);

  // Description:
  // Set the layout strategy to a tree layout. Radial indicates whether to
  // do a radial or standard top-down tree layout. The angle parameter is the
  // angular distance spanned by the tree. Leaf spacing is a
  // value from 0 to 1 indicating how much of the radial layout should be
  // allocated to leaf nodes (as opposed to between tree branches). The log spacing value is a
  // non-negative value where > 1 will create expanding levels, < 1 will create
  // contracting levels, and = 1 makes all levels the same size. See
  // vtkTreeLayoutStrategy for more information.
  virtual void SetLayoutStrategyToTree(
    bool radial = false,
    double angle = 90,
    double leafSpacing = 0.9,
    double logSpacing = 1.0);

  // Description:
  // Set the layout strategy to a cosmic tree layout. nodeSizeArrayName is
  // the array used to size the circles (default is NULL, which makes leaf
  // nodes the same size). sizeLeafNodesOnly only uses the leaf node sizes,
  // and computes the parent size as the sum of the child sizes (default true).
  // layoutDepth stops layout at a certain depth (default is 0, which does the
  // entire tree). layoutRoot is the vertex that will be considered the root
  // node of the layout (default is -1, which will use the tree's root).
  // See vtkCosmicTreeLayoutStrategy for more information.
  virtual void SetLayoutStrategyToCosmicTree(
    const char* nodeSizeArrayName = 0,
    bool sizeLeafNodesOnly = true,
    int layoutDepth = 0,
    vtkIdType layoutRoot = -1);

  // Description:
  // Set/get the graph layout strategy.
  virtual void SetEdgeLayoutStrategy(vtkEdgeLayoutStrategy* strategy);
  virtual vtkEdgeLayoutStrategy* GetEdgeLayoutStrategy();
  void SetEdgeLayoutStrategyToArcParallel()
    { this->SetEdgeLayoutStrategy("Arc Parallel"); }
  void SetEdgeLayoutStrategyToPassThrough()
    { this->SetEdgeLayoutStrategy("Pass Through"); }

  // Description:
  // Set the edge layout strategy to a geospatial arced strategy
  // appropriate for vtkGeoView.
  virtual void SetEdgeLayoutStrategyToGeo(double explodeFactor = 0.2);

  // Description:
  // Set the edge layout strategy by name.
  virtual void SetEdgeLayoutStrategy(const char* name);
  vtkGetStringMacro(EdgeLayoutStrategyName);

  // Description:
  // Apply a theme to this representation.
  virtual void ApplyViewTheme(vtkViewTheme* theme);

  // Description:
  // Set the graph vertex glyph type.
  virtual void SetGlyphType(int type);
  virtual int GetGlyphType();

  // Description:
  // Set whether to scale vertex glyphs.
  virtual void SetScaling(bool b);
  virtual bool GetScaling();
  vtkBooleanMacro(Scaling, bool);

  // Description:
  // Set the glyph scaling array name.
  virtual void SetScalingArrayName(const char* name);
  virtual const char* GetScalingArrayName();

  // Description:
  // Vertex/edge scalar bar visibility.
  virtual void SetVertexScalarBarVisibility(bool b);
  virtual bool GetVertexScalarBarVisibility();
  virtual void SetEdgeScalarBarVisibility(bool b);
  virtual bool GetEdgeScalarBarVisibility();

  // Description:
  // Whether the current graph layout is complete.
  virtual bool IsLayoutComplete();

  // Description:
  // Performs another iteration on the graph layout.
  virtual void UpdateLayout();

protected:
  vtkRenderedGraphRepresentation();
  ~vtkRenderedGraphRepresentation();

  // Description:
  // Called by the view to add/remove this representation.
  virtual bool AddToView(vtkView* view);
  virtual bool RemoveFromView(vtkView* view);

  virtual void PrepareForRendering(vtkRenderView* view);
  
  virtual vtkSelection* ConvertSelection(vtkView* view, vtkSelection* sel);

  // Description:
  // Sets up the input connections for this representation.
  virtual void SetupInputConnections();

  //BTX
  // Description:
  // Internal filter classes.
  vtkSmartPointer<vtkApplyColors>          ApplyColors;
  vtkSmartPointer<vtkVertexDegree>         VertexDegree;
  vtkSmartPointer<vtkPolyData>             EmptyPolyData;
  vtkSmartPointer<vtkEdgeCenters>          EdgeCenters;
  vtkSmartPointer<vtkGraphToPoints>        GraphToPoints;
  vtkSmartPointer<vtkArrayMap>             VertexLabels;
  vtkSmartPointer<vtkArrayMap>             EdgeLabels;
  vtkSmartPointer<vtkArrayMap>             VertexLabelPriority;
  vtkSmartPointer<vtkArrayMap>             EdgeLabelPriority;
  vtkSmartPointer<vtkTextProperty>         VertexTextProperty;
  vtkSmartPointer<vtkTextProperty>         EdgeTextProperty;
  vtkSmartPointer<vtkArrayMap>             VertexIcons;
  vtkSmartPointer<vtkArrayMap>             EdgeIcons;
  vtkSmartPointer<vtkArrayMap>             VertexIconPriority;
  vtkSmartPointer<vtkArrayMap>             EdgeIconPriority;
  vtkSmartPointer<vtkGraphLayout>          Layout;
  vtkSmartPointer<vtkPerturbCoincidentVertices> Coincident;
  vtkSmartPointer<vtkEdgeLayout>           EdgeLayout;
  vtkSmartPointer<vtkGraphToPolyData>      GraphToPoly;
  vtkSmartPointer<vtkPolyDataMapper>       EdgeMapper;
  vtkSmartPointer<vtkActor>                EdgeActor;
  vtkSmartPointer<vtkGraphToGlyphs>        VertexGlyph;
  vtkSmartPointer<vtkPolyDataMapper>       VertexMapper;
  vtkSmartPointer<vtkActor>                VertexActor;
  vtkSmartPointer<vtkGraphToGlyphs>        OutlineGlyph;
  vtkSmartPointer<vtkPolyDataMapper>       OutlineMapper;
  vtkSmartPointer<vtkActor>                OutlineActor;
  vtkSmartPointer<vtkScalarBarWidget>      VertexScalarBar;
  vtkSmartPointer<vtkScalarBarWidget>      EdgeScalarBar;
  //ETX

  vtkSetStringMacro(VertexColorArrayNameInternal);
  vtkGetStringMacro(VertexColorArrayNameInternal);
  char* VertexColorArrayNameInternal;

  vtkSetStringMacro(EdgeColorArrayNameInternal);
  vtkGetStringMacro(EdgeColorArrayNameInternal);
  char* EdgeColorArrayNameInternal;

  vtkSetStringMacro(ScalingArrayNameInternal);
  vtkGetStringMacro(ScalingArrayNameInternal);
  char* ScalingArrayNameInternal;

  vtkSetStringMacro(LayoutStrategyName);
  char* LayoutStrategyName;
  vtkSetStringMacro(EdgeLayoutStrategyName);
  char* EdgeLayoutStrategyName;

private:
  vtkRenderedGraphRepresentation(const vtkRenderedGraphRepresentation&); // Not implemented
  void operator=(const vtkRenderedGraphRepresentation&);   // Not implemented
};

#endif

