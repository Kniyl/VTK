/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkGLPolyDataMapper.cxx
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
#include "vtkGLRenderer.h"
#include "vtkPolygon.h"
#include "vtkTriangle.h"
#include "vtkPolyData.h"
#include "vtkGLPolyDataMapper.h"

// Description:
// Get the lmcolor property, this is a pretty important little 
// function.  It determines how vertex colors will be handled  
// in gl.  When a primitive has vertex colors it will use this 
// method to determine what lmcolor mode to set.               
int vtkGLPolyDataMapper::GetLmcolorMode(vtkProperty *prop)
{
  if (prop->GetAmbient() > prop->GetDiffuse())
    {
    return LMC_AMBIENT;
    }
  else
    {
    return LMC_DIFFUSE;
    }
}

//
// Receives from Actor -> maps data to primitives
//
void vtkGLPolyDataMapper::Render(vtkRenderer *ren, vtkActor *act)
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
  // if something has changed regenrate colors and display lists
  // if required
  //
  if ( this->GetMTime() > this->BuildTime || 
       input->GetMTime() > this->BuildTime || 
       this->LookupTable->GetMTime() > this->BuildTime ||
       act->GetProperty()->GetMTime() > this->BuildTime)
    {
    // sets this->Colors as side effect
    this->GetColors();
    this->BuildTime.Modified();
    }

  // want to draw the primitives here
  this->Draw(ren,act);
}

// Description:
// Load poly data into gl graphics library.
void vtkGLPolyDataMapper::Draw(vtkRenderer *vtkNotUsed(aren), vtkActor *act)
{
  int npts, idx[3], rep, j, interpolation;
  float fclr[4], polyNorm[3], tran;
  short clr[4];
  void (*bgn_func[4])(),(*end_func[4])();
  void (*aBgn_func)(),(*aEnd_func)();
  vtkProperty *prop;
  vtkPoints *p;
  vtkCellArray *prims[4], *aPrim;
  vtkScalars *c=NULL;
  vtkNormals *n;
  unsigned char *rgba;
  int *pts;
  vtkTCoords *t;
  int tDim, primType;
  vtkPolyData *input = (vtkPolyData *)this->Input;
  
  // get the property 
  prop = act->GetProperty();

  // get the transparency 
  tran = prop->GetOpacity();
  
  // if the primitives are invisable then get out of here 
  if (tran <= 0.0) return;

  // get the representation (e.g., surface / wireframe / points)
  rep = prop->GetRepresentation();

  switch (rep) 
    {
    case VTK_POINTS:
      bgn_func[0] = bgnpoint;
      bgn_func[1] = bgnpoint;
      bgn_func[2] = bgnpoint;
      bgn_func[3] = bgnpoint;
      end_func[0] = endpoint;
      end_func[1] = endpoint;
      end_func[2] = endpoint;
      end_func[3] = endpoint;
      break;
    case VTK_WIREFRAME:
      bgn_func[0] = bgnpoint;
      bgn_func[1] = bgnline;
      bgn_func[2] = bgnline;
      bgn_func[3] = bgnclosedline;
      end_func[0] = endpoint;
      end_func[1] = endline;
      end_func[2] = endline;
      end_func[3] = endclosedline;
      break;
    case VTK_SURFACE:
      bgn_func[0] = bgnpoint;
      bgn_func[1] = bgnline;
      bgn_func[2] = bgntmesh;
      bgn_func[3] = bgnpolygon;
      end_func[0] = endpoint;
      end_func[1] = endline;
      end_func[2] = endtmesh;
      end_func[3] = endpolygon;
      break;
    default: 
      vtkErrorMacro(<< "Bad glr_poly representation sent\n");
      bgn_func[0] = bgnpoint;
      bgn_func[1] = bgnline;
      bgn_func[2] = bgntmesh;
      bgn_func[3] = bgnpolygon;
      end_func[0] = endpoint;
      end_func[1] = endline;
      end_func[2] = endtmesh;
      end_func[3] = endpolygon;
      break;
    }

  // get the shading interpolation 
  interpolation = prop->GetInterpolation();

  // and draw the display list
  p = input->GetPoints();
  if ( this->Colors )
    {
    c = this->Colors;
    c->InitColorTraversal(tran, this->LookupTable, this->ColorMode);
    }
    
  prims[0] = input->GetVerts();
  prims[1] = input->GetLines();
  prims[2] = input->GetStrips();
  prims[3] = input->GetPolys();

  t = input->GetPointData()->GetTCoords();
  if ( t ) 
    {
    tDim = t->GetGetNumberOfComponents();
    if (tDim != 2)
      {
      vtkDebugMacro(<< "Currently only 2d textures are supported.\n");
      t = NULL;
      }
    }

  n = input->GetPointData()->GetNormals();
  if (interpolation == VTK_FLAT) n = 0;

  // if we are doing vertex colors then set lmcolor to adjust 
  // the current materials ambient and diffuse values using   
  // vertex color commands otherwise tell it not to.          
  if (this->Colors)
    {
    lmcolor(this->GetLmcolorMode(prop));
    }
  else 
    {
    lmcolor(LMC_NULL);
    }
  
  for (primType = 0; primType < 4; primType++)
    {
    aPrim = prims[primType];
    aBgn_func = bgn_func[primType];
    aEnd_func = end_func[primType];

    // for lines or points
    if (primType < 2 && !c)
      {
      float *bg_color;
      float ambient;
      // if a line is being drawn without normals and with the  
      // ambient intensity set to zero, then lets pretend that  
      // the ambient intensity is 1.0 because otherwise the line
      // would either not show up or be screwed up              
      ambient = prop->GetAmbient();
      if (ambient <= 0.0)
	{
	// get the color from the property and set it 
	bg_color = prop->GetColor();
	fclr[0] = bg_color[0]; 
	fclr[1] = bg_color[1]; 
	fclr[2] = bg_color[2];
	fclr[3]  = tran;
	lmcolor(LMC_COLOR);
	bgnpoint();
	c4f(fclr);
	endpoint();
	}
      else
	lmcolor(LMC_NULL);
      }
    
    for (aPrim->InitTraversal(); aPrim->GetNextCell(npts,pts); )
      { 
      (*aBgn_func)();
      
      if ((primType > 1) && (!n))
        {
        if ( primType == 3 ) vtkPolygon::ComputeNormal(p,npts,pts,polyNorm);
	else vtkTriangle::ComputeNormal(p,3,pts,polyNorm);
        }
      
      for (j = 0; j < npts; j++) 
	{
	if (c) 
	  {
	  rgba = c->GetColor(pts[j]);
	  clr[0] = rgba[0]; 
	  clr[1] = rgba[1]; 
	  clr[2] = rgba[2];
	  clr[3] = rgba[3];
	  c4s(clr);
	  }
	
	if (t)
	  {
	  t2f (t->GetTCoord(pts[j]));
	  }
	
	if (n) 
	  {
	  n3f(n->GetNormal(pts[j]));
	  }
	else 
	  {
	  if (primType == 3) 
	    {
	    n3f(polyNorm);
	    }
	  if (primType == 2)
	    {
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
	    n3f(polyNorm);
	    }
	  }
	
	v3f(p->GetPoint(pts[j]));
	}
      (*aEnd_func)();

      // if its wireframe, then draw the top and bottom edges
      // of any tstrips
      if (primType == 2 && rep == VTK_WIREFRAME) 
	{
	// draw first line
	bgnline();
	for (j = 0; j < npts; j += 2) 
	  {
	  if (c) 
	    {
	    rgba = c->GetColor(pts[j]);
	    clr[0] = rgba[0]; 
	    clr[1] = rgba[1]; 
	    clr[2] = rgba[2];
            clr[3] = rgba[3];
	    c4s(clr);
	    }
	  
	  if (n) 
	    {
	    n3f(n->GetNormal(pts[j]));
	    }
	  else 
	    {
	    if ( j && j < (npts-1) )
	      {
              idx[0] = pts[j-1]; idx[1] = pts[j]; idx[2] = pts[j+1]; 
              vtkTriangle::ComputeNormal(p, 3, idx, polyNorm);
	      }
	    n3f(polyNorm);
	    }
	  
	  if (t)
	    {
	    t2f (t->GetTCoord(pts[j]));
	    }
	  
	  v3f(p->GetPoint(pts[j]));
	  }
	endline();
	
	// draw second line
	bgnline();
	for (j = 1; j < npts; j += 2) 
	  {
	  if (c) 
	    {
	    rgba = c->GetColor(pts[j]);
	    clr[0] = rgba[0]; 
	    clr[1] = rgba[1]; 
	    clr[2] = rgba[2];
            clr[3] = rgba[3];
	    c4s(clr);
	    }
	  
	  if (n) 
	    {
	    n3f(n->GetNormal(pts[j]));
	    }
	  else 
	    {
	    if (j < npts-1)
	      {
              idx[0] = pts[j+1]; idx[1] = pts[j]; idx[2] = pts[j-1]; 
              vtkTriangle::ComputeNormal(p, 3, idx, polyNorm);
	      }
	    n3f(polyNorm);
	    }
	  if (t)
	    {
	    t2f (t->GetTCoord(pts[j]));
	    }
	  v3f(p->GetPoint(pts[j]));
	  }
	endline();
	}
      }
    }
}
