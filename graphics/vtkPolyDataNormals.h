/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkPolyDataNormals.h
  Language:  C++
  Date:      $Date$
  Version:   $Revision$


Copyright (c) 1993-2000 Ken Martin, Will Schroeder, Bill Lorensen.

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
// .NAME vtkPolyDataNormals - compute normals for polygonal mesh
// .SECTION Description
// vtkPolyDataNormals is a filter that computes point normals for a polygonal 
// mesh. The filter can reorder polygons to insure consistent orientation
// across polygon neighbors. Sharp edges can be split and points duplicated
// with separate normals to give crisp (rendered) surface definition. It is
// also possible to globally flip the normal orientation.
//<P>
// The algorithm works by determining normals for each polygon and then
// averaging them at shared points. When sharp edges are present, the edges
// are split and new points generated to prevent blurry edges (due to 
// Gouraud shading).

// .SECTION Caveats
// Normals are computed only for polygons and triangle strips. Normals are
// not computed for lines or vertices.
//
// Triangle strips are broken up into triangle polygons. You may want to 
// restrip the triangles.

#ifndef __vtkPolyDataNormals_h
#define __vtkPolyDataNormals_h

#include "vtkPolyDataToPolyDataFilter.h"

class VTK_EXPORT vtkPolyDataNormals : public vtkPolyDataToPolyDataFilter
{
public:
  vtkTypeMacro(vtkPolyDataNormals,vtkPolyDataToPolyDataFilter);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Construct with feature angle=30, splitting and consistency turned on, 
  // flipNormals turned off, and non-manifold traversal turned on.
  // ComputePointNormals is on and ComputeCellNormals is off.
  static vtkPolyDataNormals *New();

  // Description:
  // Specify the angle that defines a sharp edge. If the difference in
  // angle across neighboring polygons is greater than this value, the
  // shared edge is considered "sharp".
  vtkSetClampMacro(FeatureAngle,float,0.0,180.0);
  vtkGetMacro(FeatureAngle,float);

  // Description:
  // Turn on/off the splitting of sharp edges.
  vtkSetMacro(Splitting,int);
  vtkGetMacro(Splitting,int);
  vtkBooleanMacro(Splitting,int);

  // Description:
  // Turn on/off the enforcement of consistent polygon ordering.
  vtkSetMacro(Consistency,int);
  vtkGetMacro(Consistency,int);
  vtkBooleanMacro(Consistency,int);

  // Description:
  // Turn on/off the computation of Point Normals.
  vtkSetMacro(ComputePointNormals,int);
  vtkGetMacro(ComputePointNormals,int);
  vtkBooleanMacro(ComputePointNormals,int);

  // Description:
  // Turn on/off the computation of Cell Normals.
  vtkSetMacro(ComputeCellNormals,int);
  vtkGetMacro(ComputeCellNormals,int);
  vtkBooleanMacro(ComputeCellNormals,int);

  // Description:
  // Turn on/off the global flipping of normal orientation. Flipping
  // reverves the meaning of front and back for Frontface and Backface
  // culling in vtkProperty.  Flipping modifies both the normal
  // direction and the order of a cell's points.
  vtkSetMacro(FlipNormals,int);
  vtkGetMacro(FlipNormals,int);
  vtkBooleanMacro(FlipNormals,int);

  // Description:
  // Control the depth of recursion used in this algorithm. (Some systems
  // have limited stack depth.)
  vtkSetClampMacro(MaxRecursionDepth,int,10,VTK_LARGE_INTEGER);
  vtkGetMacro(MaxRecursionDepth,int);

  // Description:
  // Turn on/off traversal across non-manifold edges. This will prevent
  // problems where the consistency of polygonal ordering is corrupted due
  // to topological loops.
  vtkSetMacro(NonManifoldTraversal,int);
  vtkGetMacro(NonManifoldTraversal,int);
  vtkBooleanMacro(NonManifoldTraversal,int);

protected:
  vtkPolyDataNormals();
  ~vtkPolyDataNormals() {};
  vtkPolyDataNormals(const vtkPolyDataNormals&) {};
  void operator=(const vtkPolyDataNormals&) {};

  // Usual data generation method
  void Execute();

  float FeatureAngle;
  int Splitting;
  int Consistency;
  int FlipNormals;
  int MaxRecursionDepth;
  int NonManifoldTraversal;
  int ComputePointNormals;
  int ComputeCellNormals;

  int Mark;
  int NumFlips;
  // Returns 1 if neighbor cells need to be checked, 0 otherwise.
  // edgeNeighbors is temp storage for this routine to use. We do not
  // access its values outside this routine. This is to avoid
  // New/Delete's within the routine.
  int TraverseAndOrder(int cellId, vtkIdList *edgeNeighbors,
		       int *Visited,
		       vtkPolyData *OldMesh, vtkPolyData *NewMesh);

  // edgeNeighbors is temp storage for this routine to use. We do not
  // access its values outside this routine. This is to avoid
  // New/Delete's within the routine.
  void MarkAndReplace(int cellId, int n, int replacement,
		      vtkNormals *PolyNormals, vtkIdList *edgeNeighbors,
		      int *Visited, vtkIdList *Map,
		      vtkPolyData *OldMesh, vtkPolyData *NewMesh,
		      float CosAngle);
};

#endif


