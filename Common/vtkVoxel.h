/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkVoxel.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkVoxel - a cell that represents a 3D orthogonal parallelepiped
// .SECTION Description
// vtkVoxel is a concrete implementation of vtkCell to represent a 3D
// orthogonal parallelepiped. Unlike vtkHexahedron, vtkVoxel has interior
// angles of 90 degrees, and sides are parallel to coordinate axes. This
// results in large increases in computational performance.

#ifndef __vtkVoxel_h
#define __vtkVoxel_h

#include "vtkCell3D.h"

class vtkLine;
class vtkPixel;

class VTK_COMMON_EXPORT vtkVoxel : public vtkCell3D
{
public:
  static vtkVoxel *New();
  vtkTypeRevisionMacro(vtkVoxel,vtkCell3D);

  // Description:
  // See vtkCell3D API for description of these methods.
  virtual void GetEdgePoints(int edgeId, int* &pts);
  virtual void GetFacePoints(int faceId, int* &pts);
  virtual double *GetParametricCoords();

  // Description:
  // See the vtkCell API for descriptions of these methods.
  int GetCellType() {return VTK_VOXEL;}
  int GetCellDimension() {return 3;}
  int GetNumberOfEdges() {return 12;}
  int GetNumberOfFaces() {return 6;}
  vtkCell *GetEdge(int edgeId);
  vtkCell *GetFace(int faceId);
  int CellBoundary(int subId, double pcoords[3], vtkIdList *pts);
  void Contour(double value, vtkDataArray *cellScalars, 
               vtkPointLocator *locator, vtkCellArray *verts, 
               vtkCellArray *lines, vtkCellArray *polys,
               vtkPointData *inPd, vtkPointData *outPd,
               vtkCellData *inCd, vtkIdType cellId, vtkCellData *outCd);
  int EvaluatePosition(double x[3], double* closestPoint,
                       int& subId, double pcoords[3],
                       double& dist2, double *weights);
  void EvaluateLocation(int& subId, double pcoords[3], double x[3],
                        double *weights);
  int IntersectWithLine(double p1[3], double p2[3], double tol, double& t,
                        double x[3], double pcoords[3], int& subId);
  int Triangulate(int index, vtkIdList *ptIds, vtkPoints *pts);
  void Derivatives(int subId, double pcoords[3], double *values, 
                   int dim, double *derivs);

  // Description:
  // Voxel specific methods for interpolation and derivatives.
  static void InterpolationFunctions(double pcoords[3], double weights[8]);
  static void InterpolationDerivs(double pcoords[3], double derivs[24]);
  static int *GetEdgeArray(int edgeId);
  static int *GetFaceArray(int faceId);

protected:
  vtkVoxel();
  ~vtkVoxel();

  vtkLine *Line;
  vtkPixel *Pixel;
  
private:
  vtkVoxel(const vtkVoxel&);  // Not implemented.
  void operator=(const vtkVoxel&);  // Not implemented.
};

#endif


