/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkTriangle.h
  Language:  C++
  Date:      $Date$
  Version:   $Revision$


Copyright (c) 1993-1998 Ken Martin, Will Schroeder, Bill Lorensen.

This software is copyrighted by Ken Martin, Will Schroeder and Bill Lorensen.
The following terms apply to all files associated with the software unless
explicitly disclaimed in individual files. This copyright specifically does
not apply to the related textbook "The Visualization Toolkit" ISBN
013199837-4 published by Prentice Hall which is covered by its own copyright.

The authors hereby grant permission to use, copy, and distribute this
software and its documentation for any purpose, provided that existing
copyright notices are retained in all copies and that this notice is included
verbatim in any distributions. Additionally, the authors grant permission to
modify this software and its documentation for any purpose, provided that
such modifications are not distributed without the explicit consent of the
authors and that existing copyright notices are retained in all copies. Some
of the algorithms implemented by this software are patented, observe all
applicable patent law.

IN NO EVENT SHALL THE AUTHORS OR DISTRIBUTORS BE LIABLE TO ANY PARTY FOR
DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES ARISING OUT
OF THE USE OF THIS SOFTWARE, ITS DOCUMENTATION, OR ANY DERIVATIVES THEREOF,
EVEN IF THE AUTHORS HAVE BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

THE AUTHORS AND DISTRIBUTORS SPECIFICALLY DISCLAIM ANY WARRANTIES, INCLUDING,
BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
PARTICULAR PURPOSE, AND NON-INFRINGEMENT.  THIS SOFTWARE IS PROVIDED ON AN
"AS IS" BASIS, AND THE AUTHORS AND DISTRIBUTORS HAVE NO OBLIGATION TO PROVIDE
MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.


=========================================================================*/
// .NAME vtkTriangle - a cell that represents a triangle
// .SECTION Description
// vtkTriangle is a concrete implementation of vtkCell to represent a triangle
// located in 3-space.

#ifndef __vtkTriangle_h
#define __vtkTriangle_h

#include "vtkCell.h"
#include "vtkMath.h"
#include "vtkLine.h"

class VTK_EXPORT vtkTriangle : public vtkCell
{
public:
  vtkTriangle();
  ~vtkTriangle();
  static vtkTriangle *New() {return new vtkTriangle;};
  const char *GetClassName() {return "vtkTriangle";};

  // Description:
  // Create a new cell and copy this triangle's information into the
  // cell. Returns a poiner to the new cell created.
  vtkCell *MakeObject();

  // Description:
  // Get the edge specified by edgeId (range 0 to 2) and return that edge's
  // coordinates.
  vtkCell *GetEdge(int edgeId);

  // Description:
  // See the vtkCell API for descriptions of these methods.
  int GetCellType() {return VTK_TRIANGLE;};
  int GetCellDimension() {return 2;};
  int GetNumberOfEdges() {return 3;};
  int GetNumberOfFaces() {return 0;};
  vtkCell *GetFace(int) {return 0;};
  int CellBoundary(int subId, float pcoords[3], vtkIdList *pts);
  void Contour(float value, vtkScalars *cellScalars, 
               vtkPointLocator *locator, vtkCellArray *verts,
               vtkCellArray *lines, vtkCellArray *polys, 
               vtkPointData *inPd, vtkPointData *outPd,
               vtkCellData *inCd, int cellId, vtkCellData *outCd);
  int EvaluatePosition(float x[3], float closestPoint[3],
                       int& subId, float pcoords[3],
                       float& dist2, float *weights);
  void EvaluateLocation(int& subId, float pcoords[3], float x[3],
                        float *weights);
  int Triangulate(int index, vtkIdList *ptIds, vtkPoints *pts);
  void Derivatives(int subId, float pcoords[3], float *values, 
                   int dim, float *derivs);

  // Description:
  // Clip this triangle using scalar value provided. Like contouring, except
  // that it cuts the triangle to produce other triangles.
  void Clip(float value, vtkScalars *cellScalars, 
            vtkPointLocator *locator, vtkCellArray *polys,
            vtkPointData *inPd, vtkPointData *outPd,
            vtkCellData *inCd, int cellId, vtkCellData *outCd, int insideOut);

  // Description:
  // Plane intersection plus in/out test on triangle. The in/out test is 
  // performed using tol as the tolerance.
  int IntersectWithLine(float p1[3], float p2[3], float tol, float& t,
                        float x[3], float pcoords[3], int& subId);

  // Description:
  // Return the center of the triangle in parametric coordinates.
  int GetParametricCenter(float pcoords[3]);

  // Description:
  // Compute the center of the triangle.
  static void TriangleCenter(float p1[3], float p2[3], float p3[3], 
                             float center[3]);

  // Description:
  // Compute the area of a triangle in 3D.
  static float TriangleArea(float p1[3], float p2[3], float p3[3]);
  
  // Description:
  // Compute the circumcenter (center[3]) and radius (method return value) of
  // a triangle defined by the three points x1, x2, and x3. (Note that the
  // coordinates are 2D. 3D points can be used but the z-component will be
  // ignored.)
  static double Circumcircle(double  p1[2], double p2[2], double p3[2], 
                            double center[2]);

  // Description:
  // Given a 2D point x[2], determine the barycentric coordinates of the point.
  // Barycentric coordinates are a natural coordinate system for simplices that
  // express a position as a linear combination of the vertices. For a 
  // triangle, there are three barycentric coordinates (because there are
  // fourthree vertices), and the sum of the coordinates must equal 1. If a 
  // point x is inside a simplex, then all three coordinates will be strictly 
  // positive.  If two coordinates are zero (so the third =1), then the 
  // point x is on a vertex. If one coordinates are zero, the point x is on an 
  // edge. In this method, you must specify the vertex coordinates x1->x3. 
  // Returns 0 if triangle is degenerate.
  static int BarycentricCoords(double x[2], double  x1[2], double x2[2], 
                               double x3[2], double bcoords[3]);
  
  
  // Description:
  // Project triangle defined in 3D to 2D coordinates. Returns 0 if
  // degenerate triangle; non-zero value otherwise. Input points are x1->x3;
  // output 2D points are v1->v3.
  static int ProjectTo2D(double x1[3], double x2[3], double x3[3],
                         double v1[2], double v2[2], double v3[2]);

  // Description:
  // Compute the triangle normal from a points list, and a list of point ids
  // that index into the points list.
  static void ComputeNormal(vtkPoints *p, int numPts, int *pts, float n[3]);

  // Description:
  // Compute the triangle normal from three points.
  static void ComputeNormal(float v1[3], float v2[3], float v3[3], float n[3]);
  
  // Description:
  // Compute the triangle normal from three points (double-precision version).
  static void ComputeNormal(double v1[3], double v2[3], double v3[3], 
                            double n[3]);
  
  // Description:
  // Given a point x, determine whether it is inside (within the
  // tolerance squared, tol2) the triangle defined by the three 
  // coordinate values p1, p2, p3. Method is via comparing dot products.
  // (Note: in current implementation the tolerance only works in the
  // neighborhood of the three vertices of the triangle.
  static int PointInTriangle(float x[3], float x1[3], float x2[3], float x3[3], 
                             float tol2);

  // Description:
  // For legacy compatibility. Do not use.
  int CellBoundary(int subId, float pcoords[3], vtkIdList &pts)
    {return this->CellBoundary(subId, pcoords, &pts);}
  int Triangulate(int index, vtkIdList &ptIds, vtkPoints &pts)
    {return this->Triangulate(index, &ptIds, &pts);}
  

protected:
  vtkLine *Line;

};

inline int vtkTriangle::GetParametricCenter(float pcoords[3])
{
  pcoords[0] = pcoords[1] = 0.333; pcoords[2] = 0.0;
  return 0;
}

inline void vtkTriangle::ComputeNormal(float v1[3], float v2[3], 
                                       float v3[3], float n[3])
{
  float length, ax, ay, az, bx, by, bz;

  // order is important!!! maintain consistency with triangle vertex order 
  ax = v3[0] - v2[0]; ay = v3[1] - v2[1]; az = v3[2] - v2[2];
  bx = v1[0] - v2[0]; by = v1[1] - v2[1]; bz = v1[2] - v2[2];

  n[0] = (ay * bz - az * by);
  n[1] = (az * bx - ax * bz);
  n[2] = (ax * by - ay * bx);

  if ( (length = sqrt((n[0]*n[0] + n[1]*n[1] + n[2]*n[2]))) != 0.0 )
    {
    n[0] /= length;
    n[1] /= length;
    n[2] /= length;
    }
}

inline void vtkTriangle::ComputeNormal(double v1[3], double v2[3], 
                                       double v3[3], double n[3])
{
  double length, ax, ay, az, bx, by, bz;

  // order is important!!! maintain consistency with triangle vertex order 
  ax = v3[0] - v2[0]; ay = v3[1] - v2[1]; az = v3[2] - v2[2];
  bx = v1[0] - v2[0]; by = v1[1] - v2[1]; bz = v1[2] - v2[2];

  n[0] = (ay * bz - az * by);
  n[1] = (az * bx - ax * bz);
  n[2] = (ax * by - ay * bx);

  if ( (length = sqrt((n[0]*n[0] + n[1]*n[1] + n[2]*n[2]))) != 0.0 )
    {
    n[0] /= length;
    n[1] /= length;
    n[2] /= length;
    }
}

inline void vtkTriangle::TriangleCenter(float p1[3], float p2[3], float p3[3],
                                       float center[3])
{
  center[0] = (p1[0]+p2[0]+p3[0]) / 3.0;
  center[1] = (p1[1]+p2[1]+p3[1]) / 3.0;
  center[2] = (p1[2]+p2[2]+p3[2]) / 3.0;
}

inline float vtkTriangle::TriangleArea(float p1[3], float p2[3], float p3[3])
{
  float a,b,c;
  a = vtkMath::Distance2BetweenPoints(p1,p2);
  b = vtkMath::Distance2BetweenPoints(p2,p3);
  c = vtkMath::Distance2BetweenPoints(p3,p1);
  return (0.25* sqrt(fabs((double)4.0*a*c - (a-b+c)*(a-b+c))));
} 

#endif


