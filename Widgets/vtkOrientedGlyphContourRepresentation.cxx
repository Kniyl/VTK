/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkOrientedGlyphContourRepresentation.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkOrientedGlyphContourRepresentation.h"
#include "vtkCleanPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkObjectFactory.h"
#include "vtkProperty.h"
#include "vtkAssemblyPath.h"
#include "vtkMath.h"
#include "vtkInteractorObserver.h"
#include "vtkLine.h"
#include "vtkCoordinate.h"
#include "vtkGlyph3D.h"
#include "vtkCursor2D.h"
#include "vtkCylinderSource.h"
#include "vtkPolyData.h"
#include "vtkPoints.h"
#include "vtkDoubleArray.h"
#include "vtkPointData.h"
#include "vtkTransformPolyDataFilter.h"
#include "vtkTransform.h"
#include "vtkCamera.h"
#include "vtkPoints.h"
#include "vtkCellArray.h"
#include "vtkFocalPlanePointPlacer.h"
#include "vtkBezierContourLineInterpolator.h"

vtkCxxRevisionMacro(vtkOrientedGlyphContourRepresentation, "1.9");
vtkStandardNewMacro(vtkOrientedGlyphContourRepresentation);

//----------------------------------------------------------------------
vtkOrientedGlyphContourRepresentation::vtkOrientedGlyphContourRepresentation()
{
  // Initialize state
  this->InteractionState = vtkContourRepresentation::Outside;

  this->CursorShape = NULL;
  this->ActiveCursorShape = NULL;

  this->HandleSize = 0.01;
  
  this->PointPlacer = vtkFocalPlanePointPlacer::New();
  this->LineInterpolator = vtkBezierContourLineInterpolator::New();
  
  // Represent the position of the cursor
  this->FocalPoint = vtkPoints::New();
  this->FocalPoint->SetNumberOfPoints(100);
  this->FocalPoint->SetNumberOfPoints(1);
  this->FocalPoint->SetPoint(0, 0.0,0.0,0.0);
  
  vtkDoubleArray *normals = vtkDoubleArray::New();
  normals->SetNumberOfComponents(3);
  normals->SetNumberOfTuples(100);
  normals->SetNumberOfTuples(1);
  double n[3] = {0,0,0};
  normals->SetTuple(0,n);
  
  // Represent the position of the cursor
  this->ActiveFocalPoint = vtkPoints::New();
  this->ActiveFocalPoint->SetNumberOfPoints(100);
  this->ActiveFocalPoint->SetNumberOfPoints(1);
  this->ActiveFocalPoint->SetPoint(0, 0.0,0.0,0.0);
  
  vtkDoubleArray *activeNormals = vtkDoubleArray::New();
  activeNormals->SetNumberOfComponents(3);
  activeNormals->SetNumberOfTuples(100);
  activeNormals->SetNumberOfTuples(1);
  activeNormals->SetTuple(0,n);
  
  this->FocalData = vtkPolyData::New();
  this->FocalData->SetPoints(this->FocalPoint);
  this->FocalData->GetPointData()->SetNormals(normals);  
  normals->Delete();

  this->ActiveFocalData = vtkPolyData::New();
  this->ActiveFocalData->SetPoints(this->ActiveFocalPoint);
  this->ActiveFocalData->GetPointData()->SetNormals(activeNormals);  
  activeNormals->Delete();
  
  this->Glypher = vtkGlyph3D::New();
  this->Glypher->SetInput(this->FocalData);
  this->Glypher->SetVectorModeToUseNormal();
  this->Glypher->OrientOn();
  this->Glypher->ScalingOn();
  this->Glypher->SetScaleModeToDataScalingOff();
  this->Glypher->SetScaleFactor(1.0);

  this->ActiveGlypher = vtkGlyph3D::New();
  this->ActiveGlypher->SetInput(this->ActiveFocalData);
  this->ActiveGlypher->SetVectorModeToUseNormal();
  this->ActiveGlypher->OrientOn();
  this->ActiveGlypher->ScalingOn();
  this->ActiveGlypher->SetScaleModeToDataScalingOff();
  this->ActiveGlypher->SetScaleFactor(1.0);

  // The transformation of the cursor will be done via vtkGlyph3D
  // By default a vtkCursor2D will be used to define the cursor shape
  vtkCursor2D *cursor2D = vtkCursor2D::New();
  cursor2D->AllOff();
  cursor2D->PointOn();
  cursor2D->Update();
  this->SetCursorShape( cursor2D->GetOutput() );
  cursor2D->Delete();

  vtkCylinderSource *cylinder = vtkCylinderSource::New();
  cylinder->SetResolution(64);
  cylinder->SetRadius(0.5);
  cylinder->SetHeight(0.0);
  cylinder->CappingOff();
  cylinder->SetCenter(0,0,0);

  vtkCleanPolyData* clean = vtkCleanPolyData::New();
  clean->PointMergingOn();
  clean->CreateDefaultLocator();
  clean->SetInputConnection(0,cylinder->GetOutputPort(0));

  vtkTransform *t = vtkTransform::New();
  t->RotateZ(90.0);

  vtkTransformPolyDataFilter *tpd = vtkTransformPolyDataFilter::New();
  tpd->SetInputConnection( 0, clean->GetOutputPort(0) );
  tpd->SetTransform( t );
  clean->Delete();
  cylinder->Delete();
  
  tpd->Update();
  this->SetActiveCursorShape(tpd->GetOutput());
  tpd->Delete();
  t->Delete();
  
  this->Glypher->SetSource(this->CursorShape);
  this->ActiveGlypher->SetSource(this->ActiveCursorShape);

  this->Mapper = vtkPolyDataMapper::New();
  this->Mapper->SetInput(this->Glypher->GetOutput());
  this->Mapper->SetResolveCoincidentTopologyToPolygonOffset();
  this->Mapper->ScalarVisibilityOff();
  this->Mapper->ImmediateModeRenderingOn();

  this->ActiveMapper = vtkPolyDataMapper::New();
  this->ActiveMapper->SetInput(this->ActiveGlypher->GetOutput());
  this->ActiveMapper->SetResolveCoincidentTopologyToPolygonOffset();
  this->ActiveMapper->ScalarVisibilityOff();
  this->ActiveMapper->ImmediateModeRenderingOn();

  // Set up the initial properties
  this->CreateDefaultProperties();

  this->Actor = vtkActor::New();
  this->Actor->SetMapper(this->Mapper);
  this->Actor->SetProperty(this->Property);

  this->ActiveActor = vtkActor::New();
  this->ActiveActor->SetMapper(this->ActiveMapper);
  this->ActiveActor->SetProperty(this->ActiveProperty);

  this->Lines = vtkPolyData::New();
  this->LinesMapper = vtkPolyDataMapper::New();
  this->LinesMapper->SetInput( this->Lines );
  
  this->LinesActor = vtkActor::New();
  this->LinesActor->SetMapper( this->LinesMapper );
  this->LinesActor->SetProperty( this->LinesProperty );
  
  this->InteractionOffset[0] = 0.0;
  this->InteractionOffset[1] = 0.0;
}

//----------------------------------------------------------------------
vtkOrientedGlyphContourRepresentation::~vtkOrientedGlyphContourRepresentation()
{
  this->FocalPoint->Delete();
  this->FocalData->Delete();

  this->ActiveFocalPoint->Delete();
  this->ActiveFocalData->Delete();
  
  this->SetCursorShape( NULL );
  this->SetActiveCursorShape( NULL );

  this->Glypher->Delete();
  this->Mapper->Delete();
  this->Actor->Delete();

  this->ActiveGlypher->Delete();
  this->ActiveMapper->Delete();
  this->ActiveActor->Delete();
  
  this->Lines->Delete();
  this->LinesMapper->Delete();
  this->LinesActor->Delete();
  
  this->Property->Delete();
  this->ActiveProperty->Delete();
  this->LinesProperty->Delete();
}

//----------------------------------------------------------------------
void vtkOrientedGlyphContourRepresentation::SetCursorShape(vtkPolyData *shape)
{
  if ( shape != this->CursorShape )
    {
    if ( this->CursorShape )
      {
      this->CursorShape->Delete();
      }
    this->CursorShape = shape;
    if ( this->CursorShape )
      {
      this->CursorShape->Register(this);
      }
    if ( this->CursorShape )
      {
      this->Glypher->SetSource(this->CursorShape);
      }
    this->Modified();
    }
}

//----------------------------------------------------------------------
vtkPolyData *vtkOrientedGlyphContourRepresentation::GetCursorShape()
{
  return this->CursorShape;
}

//----------------------------------------------------------------------
void vtkOrientedGlyphContourRepresentation::SetActiveCursorShape(vtkPolyData *shape)
{
  if ( shape != this->ActiveCursorShape )
    {
    if ( this->ActiveCursorShape )
      {
      this->ActiveCursorShape->Delete();
      }
    this->ActiveCursorShape = shape;
    if ( this->ActiveCursorShape )
      {
      this->ActiveCursorShape->Register(this);
      }
    if ( this->ActiveCursorShape )
      {
      this->ActiveGlypher->SetSource(this->ActiveCursorShape);
      }
    this->Modified();
    }
}

//----------------------------------------------------------------------
vtkPolyData *vtkOrientedGlyphContourRepresentation::GetActiveCursorShape()
{
  return this->ActiveCursorShape;
}

//----------------------------------------------------------------------
void vtkOrientedGlyphContourRepresentation::SetRenderer(vtkRenderer *ren)
{
  //  this->WorldPosition->SetViewport(ren);
  this->Superclass::SetRenderer(ren);
}

//-------------------------------------------------------------------------
int vtkOrientedGlyphContourRepresentation::ComputeInteractionState(int X, int Y, int vtkNotUsed(modified))
{
  
  double pos[4], xyz[3];
  this->FocalPoint->GetPoint(0,pos);
  pos[3] = 1.0;
  this->Renderer->SetWorldPoint(pos);
  this->Renderer->WorldToDisplay();
  this->Renderer->GetDisplayPoint(pos);
  
  xyz[0] = static_cast<double>(X);
  xyz[1] = static_cast<double>(Y);
  xyz[2] = pos[2];
  
  this->VisibilityOn();
  double tol2 = this->PixelTolerance * this->PixelTolerance;
  if ( vtkMath::Distance2BetweenPoints(xyz,pos) <= tol2 )
    {
    this->InteractionState = vtkContourRepresentation::Nearby;
    if ( !this->ActiveCursorShape )
      {
      this->VisibilityOff();
      }
    }
  else
    {
    this->InteractionState = vtkContourRepresentation::Outside;
    if ( !this->CursorShape )
      {
      this->VisibilityOff();
      }
    }

  return this->InteractionState;
}

//----------------------------------------------------------------------
// Record the current event position, and the rectilinear wipe position.
void vtkOrientedGlyphContourRepresentation::StartWidgetInteraction(double startEventPos[2])
{
  this->StartEventPosition[0] = startEventPos[0];
  this->StartEventPosition[1] = startEventPos[1];
  this->StartEventPosition[2] = 0.0;

  this->LastEventPosition[0] = startEventPos[0];
  this->LastEventPosition[1] = startEventPos[1];
  
  // How far is this in pixels from the position of this widget?
  // Maintain this during interaction such as translating (don't
  // force center of widget to snap to mouse position)
  
  // convert position to display coordinates
  double pos[2];
  this->GetNthNodeDisplayPosition(this->ActiveNode, pos);

  this->InteractionOffset[0] = pos[0] - startEventPos[0];
  this->InteractionOffset[1] = pos[1] - startEventPos[1];
  
}


//----------------------------------------------------------------------
// Based on the displacement vector (computed in display coordinates) and
// the cursor state (which corresponds to which part of the widget has been
// selected), the widget points are modified.
// First construct a local coordinate system based on the display coordinates
// of the widget.
void vtkOrientedGlyphContourRepresentation::WidgetInteraction(double eventPos[2])
{
  // Process the motion
  if ( this->CurrentOperation == vtkContourRepresentation::Translate )
    {
    this->Translate(eventPos);
    }

  // Book keeping
  this->LastEventPosition[0] = eventPos[0];
  this->LastEventPosition[1] = eventPos[1];
}

//----------------------------------------------------------------------
// Translate everything
void vtkOrientedGlyphContourRepresentation::Translate(double eventPos[2])
{
  double ref[3];
  
  if ( !this->GetActiveNodeWorldPosition( ref ) )
    {
    return;
    }
  
  double displayPos[2];
  displayPos[0] = eventPos[0] + this->InteractionOffset[0];
  displayPos[1] = eventPos[1] + this->InteractionOffset[1];
  
  double worldPos[3];
  double worldOrient[9];
  if ( this->PointPlacer->ComputeWorldPosition(this->Renderer,
                                               displayPos, ref, worldPos,
                                               worldOrient ) )
    {
    this->SetActiveNodeToWorldPosition(worldPos, worldOrient);
    }
  else
    {
    // I really want to track the closest point here,
    // but I am postponing this at the moment....
    }
}


//----------------------------------------------------------------------
void vtkOrientedGlyphContourRepresentation::Scale(double eventPos[2])
{
  // Get the current scale factor
  double sf = this->Glypher->GetScaleFactor();

  // Compute the scale factor
  int *size = this->Renderer->GetSize();
  double dPos = static_cast<double>(eventPos[1]-this->LastEventPosition[1]);
  sf *= (1.0 + 2.0*(dPos / size[1])); //scale factor of 2.0 is arbitrary
  
  // Scale the handle
  this->Glypher->SetScaleFactor(sf);
}

//----------------------------------------------------------------------
void vtkOrientedGlyphContourRepresentation::CreateDefaultProperties()
{
  this->Property = vtkProperty::New();
  this->Property->SetColor(1.0,1.0,1.0);
  this->Property->SetLineWidth(0.5);
  this->Property->SetPointSize(3);

  this->ActiveProperty = vtkProperty::New();
  this->ActiveProperty->SetColor(0.0,1.0,0.0);
  this->ActiveProperty->SetRepresentationToWireframe();
  this->ActiveProperty->SetAmbient(1.0);
  this->ActiveProperty->SetDiffuse(0.0);
  this->ActiveProperty->SetSpecular(0.0);
  this->ActiveProperty->SetLineWidth(1.0);
  
  this->LinesProperty = vtkProperty::New();
  this->LinesProperty->SetAmbient(1.0);
  this->LinesProperty->SetDiffuse(0.0);
  this->LinesProperty->SetSpecular(0.0);
  this->LinesProperty->SetColor(1,1,1);
  this->LinesProperty->SetLineWidth(1);
}

//----------------------------------------------------------------------
void vtkOrientedGlyphContourRepresentation::BuildLines()
{
  vtkPoints *points = vtkPoints::New();
  vtkCellArray *lines = vtkCellArray::New();
  
  int i, j;
  vtkIdType index = 0;
  
  int count = this->GetNumberOfNodes();
  for ( i = 0; i < this->GetNumberOfNodes(); i++ )
    {
    count += this->GetNumberOfIntermediatePoints(i);
    }
  
  points->SetNumberOfPoints(count);
  vtkIdType numLines;
  
  if ( this->ClosedLoop && count > 0 )
    {
    numLines = count+1;
    }
  else
    {
    numLines = count;
    }

  if ( numLines > 0 )
    {
    vtkIdType *lineIndices = new vtkIdType[numLines];
    
    double pos[3];
    for ( i = 0; i < this->GetNumberOfNodes(); i++ )
      {
      // Add the node
      this->GetNthNodeWorldPosition( i, pos );
      points->InsertPoint( index, pos );
      lineIndices[index] = index;
      index++;
      
      int numIntermediatePoints = this->GetNumberOfIntermediatePoints(i);
      
      for ( j = 0; j < numIntermediatePoints; j++ )
        {
        this->GetIntermediatePointWorldPosition( i, j, pos );
        points->InsertPoint( index, pos );
        lineIndices[index] = index;
        index++;
        }
      }
    
    if ( this->ClosedLoop )
      {
      lineIndices[index] = 0;
      }
    
    lines->InsertNextCell( numLines, lineIndices );
    delete [] lineIndices;
    }
  
  this->Lines->SetPoints( points );
  this->Lines->SetLines( lines );
  
  points->Delete();
  lines->Delete();
}

//----------------------------------------------------------------------
vtkPolyData * 
vtkOrientedGlyphContourRepresentation::GetContourRepresentationAsPolyData()
{
  // Get the points in this contour as a vtkPolyData. 
  return this->Lines; 
}

//----------------------------------------------------------------------
void vtkOrientedGlyphContourRepresentation::BuildRepresentation()
{
  // Make sure we are up to date with any changes made in the placer
  this->UpdateContour();
  
  double p1[4], p2[4];
  this->Renderer->GetActiveCamera()->GetFocalPoint(p1);
  p1[3] = 1.0;
  this->Renderer->SetWorldPoint(p1);
  this->Renderer->WorldToView();
  this->Renderer->GetViewPoint(p1);
  
  double depth = p1[2];
  double aspect[2];
  this->Renderer->ComputeAspect();
  this->Renderer->GetAspect(aspect);
  
  p1[0] = -aspect[0];
  p1[1] = -aspect[1];
  this->Renderer->SetViewPoint(p1);
  this->Renderer->ViewToWorld();
  this->Renderer->GetWorldPoint(p1);
  
  p2[0] = aspect[0];
  p2[1] = aspect[1];
  p2[2] = depth;
  p2[3] = 1.0;
  this->Renderer->SetViewPoint(p2);
  this->Renderer->ViewToWorld();
  this->Renderer->GetWorldPoint(p2);
  
  double distance = 
    sqrt( vtkMath::Distance2BetweenPoints(p1,p2) );

  int *size = this->Renderer->GetRenderWindow()->GetSize();
  double viewport[4];
  this->Renderer->GetViewport(viewport);
  
  double x, y, scale;
  
  x = size[0] * (viewport[2]-viewport[0]);
  y = size[1] * (viewport[3]-viewport[1]);
  
  scale = sqrt( x*x + y*y );
  
  
  distance = 1000* distance / scale;
  
  this->Glypher->SetScaleFactor( distance * this->HandleSize );
  this->ActiveGlypher->SetScaleFactor( distance * this->HandleSize );
    
  int numPoints = this->GetNumberOfNodes();
  
  if ( this->ActiveNode >= 0 &&
       this->ActiveNode < this->GetNumberOfNodes() )
    {
    this->FocalPoint->SetNumberOfPoints(numPoints-1);  
    this->FocalData->GetPointData()->GetNormals()->SetNumberOfTuples(numPoints-1);
    }
  else
    {
    this->FocalPoint->SetNumberOfPoints(numPoints);  
    this->FocalData->GetPointData()->GetNormals()->SetNumberOfTuples(numPoints);
    }

  int i;
  int idx = 0;
  for ( i = 0; i < numPoints; i++ )
    {
    if ( i != this->ActiveNode )  
      {
      double worldPos[3];
      double worldOrient[9];
      this->GetNthNodeWorldPosition( i, worldPos );
      this->GetNthNodeWorldOrientation( i, worldOrient );
      this->FocalPoint->SetPoint(idx, worldPos );
      this->FocalData->GetPointData()->GetNormals()->SetTuple(idx,worldOrient+6);
      idx++;
      }
    }
  
  this->FocalPoint->Modified();
  this->FocalData->GetPointData()->GetNormals()->Modified();
  this->FocalData->Modified();

  if ( this->ActiveNode >= 0 &&
       this->ActiveNode < this->GetNumberOfNodes() )
    {
      double worldPos[3];
      double worldOrient[9];
      this->GetNthNodeWorldPosition( this->ActiveNode, worldPos );
      this->GetNthNodeWorldOrientation( this->ActiveNode, worldOrient );
      this->ActiveFocalPoint->SetPoint(0, worldPos );
      this->ActiveFocalData->GetPointData()->GetNormals()->SetTuple(0,worldOrient+6);

      this->ActiveFocalPoint->Modified();
      this->ActiveFocalData->GetPointData()->GetNormals()->Modified();
      this->ActiveFocalData->Modified();
      this->ActiveActor->VisibilityOn();
    }
  else
    {
      this->ActiveActor->VisibilityOff();
    }

}

//----------------------------------------------------------------------
void vtkOrientedGlyphContourRepresentation
::BuildRepresentationFromUserSuppliedPolydata( vtkPolyData * pd )
{
  vtkPoints *points   = pd->GetPoints();
  vtkIdType nPoints = points->GetNumberOfPoints();
  if (nPoints <= 0)
    {
    return; // Yeah right.. build from nothing !
    }

  // Clear all existing nodes.
  for(unsigned int i=0;i<this->Internal->Nodes.size();i++)
    {
    for (unsigned int j=0;j<this->Internal->Nodes[i]->Points.size();j++)
      {
      delete this->Internal->Nodes[i]->Points[j];
      }
    this->Internal->Nodes[i]->Points.clear();
    delete this->Internal->Nodes[i];
    }
  this->Internal->Nodes.clear();

  vtkIdList *pointIds = pd->GetCell(0)->GetPointIds();

  // Get the worldOrient from the point placer
  double ref[3], displayPos[2], worldPos[3], worldOrient[9];
  ref[0] = 0.0; ref[1] = 0.0; ref[2] = 0.0;
  displayPos[0] = 0.0; displayPos[1] = 0.0;
  this->PointPlacer->ComputeWorldPosition(this->Renderer,
                                 displayPos, ref, worldPos, worldOrient );

  // Add nodes

  for ( vtkIdType i=0; i < nPoints; i++ )
    {
    double *p = points->GetPoint( i );
    this->AddNodeAtWorldPosition( p, worldOrient );
    }

  if ( pointIds->GetNumberOfIds() > nPoints )
    {
    this->ClosedLoopOn();
    }

  // Update the contour representation from the nodes using the line interpolator
  this->BuildRepresentation();

  // Show the contour.
  this->VisibilityOn();
}


//----------------------------------------------------------------------
void vtkOrientedGlyphContourRepresentation::GetActors(vtkPropCollection *pc)
{
  this->Actor->GetActors(pc);
  this->ActiveActor->GetActors(pc);
}

//----------------------------------------------------------------------
void vtkOrientedGlyphContourRepresentation::ReleaseGraphicsResources(vtkWindow *win)
{
  this->Actor->ReleaseGraphicsResources(win);
  this->ActiveActor->ReleaseGraphicsResources(win);
}

//----------------------------------------------------------------------
int vtkOrientedGlyphContourRepresentation::RenderOverlay(vtkViewport *viewport)
{
  int count=0;
  count += this->LinesActor->RenderOverlay(viewport);
  if ( this->Actor->GetVisibility() )
    {
    count +=  this->Actor->RenderOverlay(viewport);
    }
  if ( this->ActiveActor->GetVisibility() )
    {
    count +=  this->ActiveActor->RenderOverlay(viewport);
    }
  return count;
}

//----------------------------------------------------------------------
int vtkOrientedGlyphContourRepresentation::RenderOpaqueGeometry(vtkViewport *viewport)
{
  // Since we know RenderOpaqueGeometry gets called first, will do the
  // build here
  this->BuildRepresentation();
  
  int count=0;
  count += this->LinesActor->RenderOpaqueGeometry(viewport);
  if ( this->Actor->GetVisibility() )
    {
    count += this->Actor->RenderOpaqueGeometry(viewport);
    }
  if ( this->ActiveActor->GetVisibility() )
    {
    count += this->ActiveActor->RenderOpaqueGeometry(viewport);
    }
  return count;
}

//----------------------------------------------------------------------
int vtkOrientedGlyphContourRepresentation::RenderTranslucentGeometry(vtkViewport *viewport)
{
  int count=0;
  count += this->LinesActor->RenderTranslucentGeometry(viewport);
  if ( this->Actor->GetVisibility() )
    {
    count += this->Actor->RenderTranslucentGeometry(viewport);
    }
  if ( this->ActiveActor->GetVisibility() )
    {
    count += this->ActiveActor->RenderTranslucentGeometry(viewport);
    }
  return count;
}


//----------------------------------------------------------------------
void vtkOrientedGlyphContourRepresentation::PrintSelf(ostream& os, vtkIndent indent)
{
  //Superclass typedef defined in vtkTypeMacro() found in vtkSetGet.h
  this->Superclass::PrintSelf(os,indent);
  
}
