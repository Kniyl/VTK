/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkPlanes.h
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
// .NAME vtkPlanes - implicit function for convex set of planes
// .SECTION Description
// vtkPlanes computes the implicit function and function gradient for a set
// of planes. The planes must define a convex space.
//
// The function value is the closest first order distance of a point to the
// convex region defined by the planes. The function gradient is the plane
// normal at the function value.  Note that the normals must point outside of
// the convex region. Thus, a negative function value means that a point is
// inside the convex region.
//
// To define the planes you must create two objects: a subclass of 
// vtkPoints (e.g., vtkPoints) and a subclass of vtkNormals (e.g., 
// vtkNormals). The points define a point on the plane, and the normals
// specify plane normals.

#ifndef __vtkPlanes_h
#define __vtkPlanes_h

#include "vtkImplicitFunction.h"

class vtkPlane;
class vtkCamera;

class VTK_EXPORT vtkPlanes : public vtkImplicitFunction
{
public:
  static vtkPlanes *New() {return new vtkPlanes;};
  const char *GetClassName() {return "vtkPlanes";};
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description
  // Evaluate plane equations. Return smallest absolute value.
  float EvaluateFunction(float x[3]);

  // Description
  // Evaluate planes gradient.
  void EvaluateGradient(float x[3], float n[3]);

  // Description:
  // Specify a list of points defining points through which the planes pass.
  vtkSetObjectMacro(Points,vtkPoints);
  vtkGetObjectMacro(Points,vtkPoints);
  
  // Description:
  // Specify a list of normal vectors for the planes. There is a one-to-one
  // correspondence between plane points and plane normals.
  vtkSetObjectMacro(Normals,vtkNormals);
  vtkGetObjectMacro(Normals,vtkNormals);

  // Description:
  // Specify the planes via a camera frustrum definition.
  void SetFrustumPlanes(vtkCamera *camera);

protected:
  vtkPlanes();
  ~vtkPlanes();
  vtkPoints *Points;
  vtkNormals *Normals;
  vtkPlane *Plane;

private:
  float Planes[24];

};

#endif


