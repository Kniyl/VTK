/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkXGLPolyDataMapper.cxx
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
#include <stdlib.h>
#include <math.h>
#include "vtkXGLRenderWindow.h"
#include "vtkXGLRenderer.h"
#include "vtkXGLPolyDataMapper.h"
#include "vtkPolygon.h"
#include "vtkTriangle.h"

// Description:
// Construct empty object.
vtkXGLPolyDataMapper::vtkXGLPolyDataMapper()
{
  this->Context = NULL;
  this->PL = NULL;
  this->PL2 = NULL;
}

vtkXGLPolyDataMapper::~vtkXGLPolyDataMapper()
{
  int i;
  float *fTemp;
  
  // free old memory 
  switch (this->FL.facet_type)
    {
    case XGL_FACET_COLOR:
      if (this->FL.facets.color_facets)
	{
	delete [] this->FL.facets.color_facets;
	this->FL.facets.color_facets = NULL;
	}
      break;
    case XGL_FACET_NORMAL:
      if (this->FL.facets.normal_facets)
	{
	delete [] this->FL.facets.normal_facets;
	this->FL.facets.normal_facets = NULL;
	}
      break;
    case XGL_FACET_COLOR_NORMAL:
      if (this->FL.facets.color_normal_facets)
	{
	delete [] this->FL.facets.color_normal_facets;
	this->FL.facets.color_normal_facets = NULL;
	}
      break;
    default: break;
    }
  this->FL.facet_type = XGL_FACET_NONE;
  this->FL.num_facets = 0;
  
  if (this->PL)
    {
    for (i = 0; i < (this->NumPolys + this->NumStrips); i++)
      {
      fTemp = (float *)(this->PL[i].pts.data_f3d);
      delete [] fTemp;
      }
    delete [] this->PL; 
    this->PL = NULL;
    }

  if (this->PL2)
    {
    for (i = 0; i < (this->NumLines + this->NumVerts); i++)
      {
      fTemp = (float *)(this->PL2[i].pts.data_f3d);
      delete [] fTemp;
      }
    delete [] this->PL2; 
    this->PL2 = NULL;
    }
}


//
// Receives from Actor -> maps data to primitives
//
void vtkXGLPolyDataMapper::Render(vtkRenderer *ren, vtkActor *act)
{
  int numPts;
  vtkPolyData *input= (vtkPolyData *)this->Input;
//
// make sure that we've been properly initialized
//
  if ( input == NULL ) 
    {
    vtkErrorMacro(<< "No input!");
    return;
    }
  else
    {
    input->Update();
    numPts = input->GetNumberOfPoints();
    } 

  if (numPts == 0)
    {
    vtkDebugMacro(<< "No points!");
    return;
    }
  
  if ( this->LookupTable == NULL ) this->CreateDefaultLookupTable();

  //
  // if something has changed regenerate colors and display lists
  // if required
  //
  if ( this->GetMTime() > this->BuildTime || 
       input->GetMTime() > this->BuildTime || 
       this->LookupTable->GetMTime() > this->BuildTime)
    {
    // sets this->Colors as side effect
    this->GetColors();
    this->Build(input, this->Colors);
    this->BuildTime.Modified();
    }
   
  // if we are in immediate mode rendering we always
  // want to draw the primitives here
  this->Draw(ren,act);
}

float *vtkXGLPolyDataMapper::AddVertex(int npts, int pointSize, int *pts,
				   vtkPoints *p, vtkScalars *c, vtkTCoords *t)
{
  float *fTemp;
  int j;
  float *pPtr;
  int adder = 0;
  unsigned char *rgb;
  
  // allocate memory
  fTemp = new float [npts*pointSize];
  
  if (fTemp == NULL) 
    {
    vtkErrorMacro(<< "XglrPoly out of memory.\n");
    return NULL;
    }
  
  for (j = 0; j < npts; j++)
    {
    pPtr = p->GetPoint(pts[j]);
    adder = j*pointSize;
    fTemp[adder] = pPtr[0];
    fTemp[adder + 1] = pPtr[1];
    fTemp[adder + 2] = pPtr[2];
    adder += 3;
    
    if (c)
      {
      rgb = c->GetColor(pts[j]);
      fTemp[adder] = rgb[0]/255.0;
      fTemp[adder + 1] = rgb[1]/255.0;
      fTemp[adder + 2] = rgb[2]/255.0;
      adder += 3;
      }
    
    if (t)
      {
      pPtr = t->GetTCoord(pts[j]);
      fTemp[adder] = pPtr[0];
      fTemp[adder + 1] = pPtr[1];
      }
    }

  return fTemp;
}


float *vtkXGLPolyDataMapper::AddVertexComputeNormal(int npts, int pointSize, 
						 int *pts, vtkPoints *p, 
						 vtkScalars *c, vtkTCoords *t, 
						 float *polyNorm)
{
  float *fTemp;
  int j, idx[3];
  float *pPtr;
  int adder = 0;
  unsigned char *rgb;
  
  // allocate memory
  fTemp = new float [npts*pointSize];
  
  if (fTemp == NULL) 
    {
    vtkErrorMacro(<< "XglrPoly out of memory.\n");
    return NULL;
    }
  
  for (j = 0; j < npts; j++)
    {
    pPtr = p->GetPoint(pts[j]);
    adder = j*pointSize;
    fTemp[adder] = pPtr[0];
    fTemp[adder + 1] = pPtr[1];
    fTemp[adder + 2] = pPtr[2];
    adder += 3;
    
    if (c)
      {
      rgb = c->GetColor(pts[j]);
      fTemp[adder] = rgb[0]/255.0;
      fTemp[adder + 1] = rgb[1]/255.0;
      fTemp[adder + 2] = rgb[2]/255.0;
      adder += 3;
      }
    
    if ( j > 2)
      {
      if (j % 2)
	{
	idx[0] = pts[j-2]; idx[1] = pts[j]; idx[2] = pts[j-1]; 
	vtkTriangle::ComputeNormal(p, 3, idx, polyNorm);
	}
      else
	{
	idx[0] = pts[j-2]; idx[1] = pts[j-1]; idx[2] = pts[j]; 
	vtkTriangle::ComputeNormal(p, 3, idx, polyNorm);
	}
      }
    else if ( j == 0 )
      {
      vtkTriangle::ComputeNormal(p, 3, pts, polyNorm);
      }
    pPtr = polyNorm;

    fTemp[adder] = pPtr[0];
    fTemp[adder + 1] = pPtr[1];
    fTemp[adder + 2] = pPtr[2];
    adder += 3;
    
    if (t)
      {
      pPtr = t->GetTCoord(pts[j]);
      fTemp[adder] = pPtr[0];
      fTemp[adder + 1] = pPtr[1];
      }
    }

  return fTemp;
}

float *vtkXGLPolyDataMapper::AddVertexWithNormal(int npts, int pointSize, 
					      int *pts, vtkPoints *p, 
					      vtkScalars *c, vtkTCoords *t, 
					      vtkNormals *n, float *polyNorm)
{
  float *fTemp;
  int j;
  float *pPtr;
  int adder = 0;
  unsigned char *rgb;
  
  // allocate memory
  fTemp = new float [npts*pointSize];
  
  if (fTemp == NULL) 
    {
    vtkErrorMacro(<< "XglrPoly out of memory.\n");
    return NULL;
    }
  
  for (j = 0; j < npts; j++)
    {
    pPtr = p->GetPoint(pts[j]);
    adder = j*pointSize;
    fTemp[adder] = pPtr[0];
    fTemp[adder + 1] = pPtr[1];
    fTemp[adder + 2] = pPtr[2];
    adder += 3;
    
    if (c)
      {
      rgb = c->GetColor(pts[j]);
      fTemp[adder] = rgb[0]/255.0;
      fTemp[adder + 1] = rgb[1]/255.0;
      fTemp[adder + 2] = rgb[2]/255.0;
      adder += 3;
      }
    
    if (n) 
      {
      pPtr = n->GetNormal(pts[j]);
      }
    else 
      {
      pPtr = polyNorm;
      }
    fTemp[adder] = pPtr[0];
    fTemp[adder + 1] = pPtr[1];
    fTemp[adder + 2] = pPtr[2];
    adder += 3;
    
    if (t)
      {
      pPtr = t->GetTCoord(pts[j]);
      fTemp[adder] = pPtr[0];
      fTemp[adder + 1] = pPtr[1];
      }
    }

  return fTemp;
}

// Description:
// Build the data structure for the XGL PolyDataMapper.
void vtkXGLPolyDataMapper::Build(vtkPolyData *data, vtkScalars *c)
{
  vtkCellArray *polys, *strips, *lines, *verts;
  vtkNormals *n;
  vtkPoints *p;
  vtkTCoords *t;
  int npts;
  float polyNorm[3];
  int i;
  int *pts;
  int tDim, pointSize, pointSize2;
  Xgl_pt_type ptType, ptType2;
  float *fTemp;
  int numDataVals;
  int cellScalars = 0;
  int cellNormals = 0;
  unsigned char *rgb;
  vtkScalars *vc;
  
  
  // free old memory
  switch (this->FL.facet_type)
    {
    case XGL_FACET_COLOR:
      if (this->FL.facets.color_facets)
	{
	delete [] this->FL.facets.color_facets;
	this->FL.facets.color_facets = NULL;
	}
      break;
    case XGL_FACET_NORMAL:
      if (this->FL.facets.normal_facets)
	{
	delete [] this->FL.facets.normal_facets;
	this->FL.facets.normal_facets = NULL;
	}
      break;
    case XGL_FACET_COLOR_NORMAL:
      if (this->FL.facets.color_normal_facets)
	{
	delete [] this->FL.facets.color_normal_facets;
	this->FL.facets.color_normal_facets = NULL;
	}
      break;
    default: break;
    }
  this->FL.facet_type = XGL_FACET_NONE;
  this->FL.num_facets = 0;

  if (this->PL)
    {
    for (i = 0; i < (this->NumPolys + this->NumStrips); i++)
      {
      fTemp = (float *)(this->PL[i].pts.data_f3d);
      delete [] fTemp;
      }
    delete [] this->PL; 
    this->PL = NULL;
    }
  if (this->PL2)
    {
    for (i = 0; i < (this->NumLines + this->NumVerts); i++)
      {
      fTemp = (float *)(this->PL2[i].pts.data_f3d);
      delete [] fTemp;
      }
    delete [] this->PL2; 
    this->PL2 = NULL;
    }

  // get the data
  polys  = data->GetPolys();
  strips = data->GetStrips();
  lines  = data->GetLines();
  verts  = data->GetVerts();

  this->NumPolys  = polys->GetNumberOfCells();
  this->NumStrips = strips->GetNumberOfCells();
  this->NumLines  = lines->GetNumberOfCells();
  this->NumVerts  = verts->GetNumberOfCells();
  this->DataSize  = 
    this->NumPolys + this->NumStrips + this->NumLines + this->NumVerts;
  
  p = data->GetPoints();
  n = data->GetPointData()->GetNormals();
  if (!n && data->GetCellData()->GetNormals())
    {
    n = data->GetCellData()->GetNormals();
    cellNormals = 1;
    }
    
  if ( c )
    {
    c->InitColorTraversal(1.0, this->LookupTable, this->ColorMode);
    if ( this->ScalarMode == VTK_SCALAR_MODE_USE_CELL_DATA ||
    !data->GetPointData()->GetScalars() )	 
      {
      cellScalars = 1;
      }
    }

  t = data->GetPointData()->GetTCoords();
  if ( t ) 
    {
    tDim = t->GetNumberOfComponents();
    if (tDim != 2)
      {
      vtkDebugMacro(<< "Currently only 2d textures are supported.\n");
      t = NULL;
      }
    }

  // get memory and copy the data
  if ((cellScalars || cellNormals) && this->NumPolys)
    {
    this->FL.num_facets = this->NumPolys;
    if (cellScalars)
      {
      if (cellNormals)
	{
	this->FL.facets.color_normal_facets = 
	  new Xgl_color_normal_facet [this->NumPolys];
	this->FL.facet_type = XGL_FACET_COLOR_NORMAL;
	this->FacetSize = sizeof(Xgl_color_normal_facet);
	}
      else
	{
	this->FL.facets.color_facets = new Xgl_color_facet [this->NumPolys];
	this->FL.facet_type = XGL_FACET_COLOR;
	this->FacetSize = sizeof(Xgl_color_facet);
	}
      }
    else
      {
      this->FL.facets.normal_facets = new Xgl_normal_facet [this->NumPolys];
      this->FL.facet_type = XGL_FACET_NORMAL;
      this->FacetSize = sizeof(Xgl_normal_facet);
      }
    }

  if (this->NumPolys + this->NumStrips)
    {
    this->PL = new Xgl_pt_list [this->NumPolys + this->NumStrips];
    if (this->PL == (Xgl_pt_list *)NULL) 
      {
      vtkErrorMacro(<< "XglrPoly out of memory.\n");
      return;
      }
    }
  if (this->NumLines + this->NumVerts)
    {
    this->PL2 = new Xgl_pt_list [this->NumLines + this->NumVerts];
    if (this->PL2 == (Xgl_pt_list *)NULL) 
      {
      vtkErrorMacro(<< "XglrPoly out of memory.\n");
      return;
      }
    }
  
  if (cellNormals)
    {
    pointSize = 3;
    }
  else
    {
    pointSize = 6;
    }
  
  if (c && !cellScalars)
    {
    pointSize += 3;
    }
  if (t)
    {
    numDataVals = 2;
    pointSize += 2;
    if (c && !cellScalars) 
      {
      if (n && cellNormals)
	{
	ptType = XGL_PT_COLOR_DATA_F3D;
	}
      else
	{
	ptType = XGL_PT_COLOR_NORMAL_DATA_F3D;
	}
      ptType2 = XGL_PT_COLOR_DATA_F3D;
      }
    else   
      {
      if (n && cellNormals)
	{
	ptType = XGL_PT_DATA_F3D;
	}
      else
	{
	ptType = XGL_PT_NORMAL_DATA_F3D;
	}
      ptType2 = XGL_PT_DATA_F3D;
      }
    }
  else
    {
    numDataVals = 0;
    if (c && !cellScalars) 
      {
      if (n && cellNormals)
	{
	ptType = XGL_PT_COLOR_F3D;
	}
      else
	{
	ptType = XGL_PT_COLOR_NORMAL_F3D;
	}
      ptType2 = XGL_PT_COLOR_F3D;
      }
    else
      {
      if (n && cellNormals)
	{
	ptType = XGL_PT_F3D;
	}
      else
	{
	ptType = XGL_PT_NORMAL_F3D;
	}
      ptType2 = XGL_PT_F3D;
      }
    }
  if (n && cellNormals)
    {
    pointSize2 = pointSize;
    }
  else
    {
    pointSize2 = pointSize - 3;
    }
  

  // vertex colors in vc if any
  if (cellScalars)
    {
    vc = NULL;
    }
  else
    {
    vc = c;
    }
  
  // DON'T MESS WITH i past this point
  i = 0;
  if (polys)
    {
    for (polys->InitTraversal(); (polys->GetNextCell(npts,pts)); i++)
      { 
      switch (this->FL.facet_type)
	{
	case XGL_FACET_COLOR: 
	  rgb = c->GetColor(i);
	  this->FL.facets.color_facets[i].color.rgb.r = rgb[0]/255.0;
	  this->FL.facets.color_facets[i].color.rgb.g = rgb[1]/255.0;
	  this->FL.facets.color_facets[i].color.rgb.b = rgb[2]/255.0;
	  break;
	case XGL_FACET_NORMAL: 
	  memcpy(&(this->FL.facets.normal_facets[i].normal),
		 n->GetNormal(i), sizeof(float)*3);
	  break;
	case XGL_FACET_COLOR_NORMAL: 
	  rgb = c->GetColor(i);
	  this->FL.facets.color_normal_facets[i].color.rgb.r = rgb[0]/255.0;
	  this->FL.facets.color_normal_facets[i].color.rgb.g = rgb[1]/255.0;
	  this->FL.facets.color_normal_facets[i].color.rgb.b = rgb[2]/255.0;
	  memcpy(&(this->FL.facets.color_normal_facets[i].normal),
		 n->GetNormal(i), sizeof(float)*3);
	  break;
	default: break;
	}
      
      this->PL[i].num_pts = npts;
      this->PL[i].bbox = NULL;
      this->PL[i].pt_type = ptType;
      this->PL[i].num_data_values = numDataVals;
      
      if (!n) vtkPolygon::ComputeNormal(p,npts,pts,polyNorm);
      
      if (cellNormals)
	{
	this->PL[i].pts.data_f3d = 
	  (Xgl_pt_data_f3d *)this->AddVertex(npts,pointSize,pts,p,vc,t);
	}
      else
	{
	this->PL[i].pts.data_f3d = 
	  (Xgl_pt_data_f3d *)this->AddVertexWithNormal(npts,pointSize,pts,
						       p,vc,t,n,polyNorm);
	}
      }
    }

  if (strips)
    {
    for (strips->InitTraversal(); (strips->GetNextCell(npts,pts)); i++)
      { 
      this->PL[i].num_pts = npts;
      this->PL[i].bbox = NULL;
      this->PL[i].pt_type = ptType;
      this->PL[i].num_data_values = numDataVals;
      
      if (n) 
	{
	this->PL[i].pts.data_f3d = 
	  (Xgl_pt_data_f3d *)this->AddVertexWithNormal(npts,pointSize,pts,
						       p,c,t,n,polyNorm);
	}
      else
	{
	this->PL[i].pts.data_f3d = 
	  (Xgl_pt_data_f3d *)this->AddVertexComputeNormal(npts,pointSize,pts,
							  p,c,t,polyNorm);
	}
      }
    }
  
  // reset i here
  i = 0;
  if (lines)
    {
    for (lines->InitTraversal(); (lines->GetNextCell(npts,pts)); i++)
      { 
      this->PL2[i].num_pts = npts;
      this->PL2[i].bbox = NULL;
      this->PL2[i].pt_type = ptType2;
      this->PL2[i].num_data_values = numDataVals;
      
      this->PL2[i].pts.data_f3d = 
	(Xgl_pt_data_f3d *)this->AddVertex(npts,pointSize2,pts,p,c,t);
      }
    }

  if (verts)
    {
    for (verts->InitTraversal(); (verts->GetNextCell(npts,pts)); i++)
      { 
      this->PL2[i].num_pts = npts;
      this->PL2[i].bbox = NULL;
      this->PL2[i].pt_type = ptType2;
      this->PL2[i].num_data_values = numDataVals;
      
      this->PL2[i].pts.data_f3d = 
	(Xgl_pt_data_f3d *)this->AddVertex(npts,pointSize2,pts,p,c,t);
      }
    }
  
  return;
}

// Description:
// Catch vtkGeometry PolyDataMapper draw method and call actual method.
void vtkXGLPolyDataMapper::Draw(vtkRenderer *aren, vtkActor *act)
{
  vtkXGLRenderer *ren = (vtkXGLRenderer *)aren;
  int polygonsToRender, i;
  Xgl_color_facet *temp;
  
  if ((!this->PL) && (!this->PL2)) return;

  this->Context = *(ren->GetContext());

  // get the transparency; return if invisible
  if ( act->GetProperty()->GetOpacity() <= 0.0 ) return;
  
  polygonsToRender = this->NumPolys;
  temp = this->FL.facets.color_facets;
  
  while (polygonsToRender > 0)
    {
    if (polygonsToRender > 100)
      {
      xgl_multi_simple_polygon(this->Context, 0, &this->FL, NULL, 100,
			       &(this->PL[this->NumPolys-polygonsToRender]));
      }
    else
      {
      xgl_multi_simple_polygon(this->Context, 0, &this->FL,
			       NULL, polygonsToRender,
			       &(this->PL[this->NumPolys-polygonsToRender]));
      }
    polygonsToRender -= 100;
    this->FL.facets.color_facets = 
      this->FL.facets.color_facets + this->FacetSize/sizeof(Xgl_color_facet);
    }
  this->FL.facets.color_facets = temp;

  
  if (this->NumStrips)
    {
    for (i = 0; i < this->NumStrips; i++)
      {
      xgl_triangle_strip(this->Context, NULL, this->PL + i + this->NumPolys);
      }
    }

  if (this->NumLines)
    {
    xgl_multipolyline(this->Context, NULL, (unsigned long)(this->NumLines), 
		      this->PL2);
    }

  if (this->NumVerts)
    {
    for (i = 0; i < this->NumVerts; i++)
      {
      xgl_multimarker(this->Context, this->PL2 + i + this->NumLines);
      }
    }
}
