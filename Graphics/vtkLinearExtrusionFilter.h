/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkLinearExtrusionFilter.h
  Language:  C++
  Date:      $Date$
  Version:   $Revision$


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
// .NAME vtkLinearExtrusionFilter - sweep polygonal data creating a "skirt" from free edges and lines, and lines from vertices
// .SECTION Description
// vtkLinearExtrusionFilter is a modeling filter. It takes polygonal data as 
// input and generates polygonal data on output. The input dataset is swept
// according to some extrusion function and creates new polygonal primitives.
// These primitives form a "skirt" or swept surface. For example, sweeping a
// line results in a quadrilateral, and sweeping a triangle creates a "wedge".
//
// There are a number of control parameters for this filter. You can 
// control whether the sweep of a 2D object (i.e., polygon or triangle strip) 
// is capped with the generating geometry via the "Capping" ivar. Also, you
// can extrude in the direction of a user specified vector, towards a point,
// or in the direction of vertex normals (normals must be provided - use 
// vtkPolyDataNormals if necessary). The amount of extrusion is controlled by
// the "ScaleFactor" instance variable.
//
// The skirt is generated by locating certain topological features. Free 
// edges (edges of polygons or triangle strips only used by one polygon or
// triangle strips) generate surfaces. This is true also of lines or 
// polylines. Vertices generate lines.
//
// This filter can be used to create 3D fonts, 3D irregular bar charts,
// or to model 2 1/2D objects like punched plates. It also can be used to 
// create solid objects from 2D polygonal meshes.

// .SECTION Caveats
// Some polygonal objects have no free edges (e.g., sphere). When swept,
// this will result in two separate surfaces if capping is on, or no surface
// if capping is off.

// .SECTION See Also
// vtkRotationalExtrusionFilter

#ifndef __vtkLinearExtrusionFilter_h
#define __vtkLinearExtrusionFilter_h

#include "vtkPolyDataToPolyDataFilter.h"

#define VTK_VECTOR_EXTRUSION 1
#define VTK_NORMAL_EXTRUSION 2
#define VTK_POINT_EXTRUSION 3

class VTK_GRAPHICS_EXPORT vtkLinearExtrusionFilter : public vtkPolyDataToPolyDataFilter 
{
public:
  vtkTypeMacro(vtkLinearExtrusionFilter,vtkPolyDataToPolyDataFilter);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Create object with normal extrusion type, capping on, scale factor=1.0,
  // vector (0,0,1), and extrusion point (0,0,0).
  static vtkLinearExtrusionFilter *New();

  // Description:
  // Set/Get the type of extrusion.
  vtkSetClampMacro(ExtrusionType,int,VTK_VECTOR_EXTRUSION,VTK_POINT_EXTRUSION);
  vtkGetMacro(ExtrusionType,int);
  void SetExtrusionTypeToVectorExtrusion()
    {this->SetExtrusionType(VTK_VECTOR_EXTRUSION);};
  void SetExtrusionTypeToNormalExtrusion()
    {this->SetExtrusionType(VTK_NORMAL_EXTRUSION);};
  void SetExtrusionTypeToPointExtrusion()
    {this->SetExtrusionType(VTK_POINT_EXTRUSION);};

  // Description:
  // Turn on/off the capping of the skirt.
  vtkSetMacro(Capping,int);
  vtkGetMacro(Capping,int);
  vtkBooleanMacro(Capping,int);

  // Description:
  // Set/Get extrusion scale factor,
  vtkSetMacro(ScaleFactor,float);
  vtkGetMacro(ScaleFactor,float);

  // Description:
  // Set/Get extrusion vector. Only needs to be set if VectorExtrusion is
  // turned on.
  vtkSetVector3Macro(Vector,float);
  vtkGetVectorMacro(Vector,float,3);

  // Description:
  // Set/Get extrusion point. Only needs to be set if PointExtrusion is
  // turned on. This is the point towards which extrusion occurs.
  vtkSetVector3Macro(ExtrusionPoint,float);
  vtkGetVectorMacro(ExtrusionPoint,float,3);

protected:
  vtkLinearExtrusionFilter();
  ~vtkLinearExtrusionFilter() {};

  void Execute();
  int ExtrusionType;
  int Capping;
  float ScaleFactor;
  float Vector[3];
  float ExtrusionPoint[3];

  //BTX
  float *(vtkLinearExtrusionFilter::*ExtrudePoint)(float x[3], vtkIdType id, 
                                                   vtkDataArray *normals);
  float *ViaNormal(float x[3], vtkIdType id, vtkDataArray *normals);
  float *ViaVector(float x[3], vtkIdType id, vtkDataArray *normals=0);
  float *ViaPoint(float x[3], vtkIdType id, vtkDataArray *normals=0);
  //ETX
 
private:
  vtkLinearExtrusionFilter(const vtkLinearExtrusionFilter&);  // Not implemented.
  void operator=(const vtkLinearExtrusionFilter&);  // Not implemented.
};

#endif
