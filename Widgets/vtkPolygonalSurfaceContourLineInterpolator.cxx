/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkPolygonalSurfaceContourLineInterpolator.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkPolygonalSurfaceContourLineInterpolator.h"

#include "vtkObjectFactory.h"
#include "vtkContourRepresentation.h"
#include "vtkPolyData.h"
#include "vtkPoints.h"
#include "vtkPointData.h"
#include "vtkCellArray.h"
#include "vtkCell.h"
#include "vtkMath.h"
#include "vtkCellLocator.h"
#include "vtkPolygonalSurfacePointPlacer.h"
#include "vtkDijkstraGraphGeodesicPath.h"

vtkCxxRevisionMacro(vtkPolygonalSurfaceContourLineInterpolator,
                                            "1.3");
vtkStandardNewMacro(vtkPolygonalSurfaceContourLineInterpolator);

//----------------------------------------------------------------------
vtkPolygonalSurfaceContourLineInterpolator
::vtkPolygonalSurfaceContourLineInterpolator()
{
  this->LastInterpolatedVertexIds[0] = -1;
  this->LastInterpolatedVertexIds[1] = -1;
  this->DistanceOffset               = 0.0;
}

//----------------------------------------------------------------------
vtkPolygonalSurfaceContourLineInterpolator
::~vtkPolygonalSurfaceContourLineInterpolator()
{
}

//----------------------------------------------------------------------
int vtkPolygonalSurfaceContourLineInterpolator::UpdateNode( 
    vtkRenderer *, vtkContourRepresentation *,
    double * vtkNotUsed(node), int vtkNotUsed(idx) )
{
  return 0;
}

//----------------------------------------------------------------------
int vtkPolygonalSurfaceContourLineInterpolator::InterpolateLine(
                          vtkRenderer *,
                          vtkContourRepresentation *rep,
                          int idx1, int idx2 )
{
  vtkPolygonalSurfacePointPlacer *placer = 
    vtkPolygonalSurfacePointPlacer::SafeDownCast(rep->GetPointPlacer());
  if (!placer)
    {
    return 1;
    }

  double p1[3], p2[3], p[3];
  rep->GetNthNodeWorldPosition( idx1, p1 );
  rep->GetNthNodeWorldPosition( idx2, p2 );

  typedef vtkPolygonalSurfacePointPlacer::Node NodeType;
  NodeType *nodeBegin = placer->GetNodeAtWorldPosition(p1);
  NodeType *nodeEnd   = placer->GetNodeAtWorldPosition(p2);
  if (nodeBegin->PolyData != nodeEnd->PolyData)
    {
    return 1;
    }

  // Find the starting and ending point id's
  vtkIdType beginVertId, endVertId;
  vtkCell *cellBegin = nodeBegin->PolyData->GetCell(nodeBegin->CellId);
  vtkPoints *cellBeginPoints = cellBegin->GetPoints();

  double minDistance = VTK_DOUBLE_MAX;
  for (int i = 0; i < cellBegin->GetNumberOfPoints(); i++)
    {
    cellBeginPoints->GetPoint(i, p);
    double distance = vtkMath::Distance2BetweenPoints( p, p1 );
    if (distance < minDistance)
      {
      beginVertId = cellBegin->GetPointId(i);
      minDistance = distance;
      }
    }

  vtkCell *cellEnd   = nodeEnd->PolyData->GetCell(nodeEnd->CellId);
  vtkPoints *cellEndPoints   = cellEnd->GetPoints();

  minDistance = VTK_DOUBLE_MAX;
  for (int i = 0; i < cellEnd->GetNumberOfPoints(); i++)
    {
    cellEndPoints->GetPoint(i, p);
    double distance = vtkMath::Distance2BetweenPoints( p, p2 );
    if (distance < minDistance)
      {
      endVertId = cellEnd->GetPointId(i);
      minDistance = distance;
      }
    }

  
  if (this->LastInterpolatedVertexIds[0] != beginVertId || 
      this->LastInterpolatedVertexIds[1] != endVertId)
    {
    vtkDijkstraGraphGeodesicPath * geodesic 
          = vtkDijkstraGraphGeodesicPath::New();
    geodesic->SetInput( nodeBegin->PolyData );
    geodesic->SetStartVertex( endVertId );
    geodesic->SetEndVertex( beginVertId );
    geodesic->Update();

    vtkPolyData *pd = geodesic->GetOutput();

    // We assume there's only one cell of course
    vtkIdType npts, *pts = NULL; 
    pd->GetLines()->InitTraversal(); 
    pd->GetLines()->GetNextCell( npts, pts );

    // Get the vertex normals if there is a height offset. The offset at
    // each node of the graph is in the direction of the vertex normal.
    
    vtkIdList *vertexIds = geodesic->GetIdList();    
    double vertexNormal[3];
    vtkDataArray *vertexNormals = NULL;
    if (this->DistanceOffset != 0.0)
      {
      vertexNormals = nodeBegin->PolyData->GetPointData()->GetNormals();
      }

    for (int n = 0; n < npts; n++)
      {
      pd->GetPoint( pts[n], p );
      if (vertexNormals)
        {
        vertexNormals->GetTuple( vertexIds->GetId(n), vertexNormal );
        p[0] += vertexNormal[0] * this->DistanceOffset;
        p[1] += vertexNormal[1] * this->DistanceOffset;
        p[2] += vertexNormal[2] * this->DistanceOffset;
        rep->AddIntermediatePointWorldPosition( idx1, p );
        }
      else
        {
        rep->AddIntermediatePointWorldPosition( idx1, p ); 
        }
      }

      this->LastInterpolatedVertexIds[0] = beginVertId;
      this->LastInterpolatedVertexIds[1] = endVertId;

      geodesic->Delete();
      }
   
  return 1;
}

//----------------------------------------------------------------------
void vtkPolygonalSurfaceContourLineInterpolator::PrintSelf(
                              ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);  

  os << indent << "DistanceOffset: " << this->DistanceOffset << endl;
}
