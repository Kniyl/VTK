/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkXglrPolyMapper.hh
  Language:  C++
  Date:      $Date$
  Version:   $Revision$


Copyright (c) 1993-1995 Ken Martin, Will Schroeder, Bill Lorensen.

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
// .NAME vtkXglrPolyMapper - a PolyMapper for Suns XGL library
// .SECTION Description
// vtkXglrPolyMapper is a subclass of vtkPolyMapperDevice.
// vtkXglrPolyMapper is a geometric PolyMapper for Suns XGL rendering library.

#ifndef __vtkXglrPolyMapper_hh
#define __vtkXglrPolyMapper_hh

#include <stdlib.h>
#include "vtkPolyMapperDevice.hh"
#include "vtkPolyData.hh"
#include "vtkColorScalars.hh"
#include <xgl/xgl.h>

class vtkXglrRenderer;
class vtkPolygon;

class vtkXglrPolyMapper : public vtkPolyMapperDevice
{
public:
  vtkXglrPolyMapper();
  virtual ~vtkXglrPolyMapper();
  char *GetClassName() {return "vtkXglrPolyMapper";};

  void Build(vtkPolyData *data, vtkColorScalars *c);
  void Draw(vtkRenderer *ren, vtkActor *a);

protected:
  float *AddVertexComputeNormal(int npts, int pointSize, int *pts, 
				vtkPoints *p, vtkColorScalars *c, 
				vtkTCoords *t, vtkPolygon *p, float *polyNorm);
  float *AddVertexWithNormal(int npts, int pointSize, int *pts, 
			     vtkPoints *p, vtkColorScalars *c, 
			     vtkTCoords *t, vtkNormals *n, float *polyNorm);
  float *AddVertex(int npts, int pointSize, int *pts, vtkPoints *p, 
		   vtkColorScalars *c, vtkTCoords *t);
  Xgl_3d_ctx Context;
  Xgl_pt_list *PL;
  Xgl_pt_list *PL2; // no normals
  int   NumPolys;
  int   NumStrips;
  int   NumLines;
  int   NumVerts;
  int   DataSize;
  
};

#endif
