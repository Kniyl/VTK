/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkCell.cxx
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
#include "vtkCell.h"

// Construct cell.
vtkCell::vtkCell()
{
  this->Points = vtkPoints::New();
  this->PointIds = vtkIdList::New();
}  

vtkCell::~vtkCell()
{
  this->Points->Delete();
  this->PointIds->Delete();
}


//
// Instantiate cell from outside
//
void vtkCell::Initialize(int npts, int *pts, vtkPoints *p)
{
  this->PointIds->Reset();
  this->Points->Reset();

  for (int i=0; i<npts; i++)
    {
    this->PointIds->InsertId(i,pts[i]);
    this->Points->InsertPoint(i,p->GetPoint(pts[i]));
    }
}
 
void vtkCell::ShallowCopy(vtkCell *c)
{
  this->Points->ShallowCopy(c->Points);
  if ( this->PointIds )
    {
    this->PointIds->Delete();
    this->PointIds = c->PointIds;
    this->PointIds->Register(this);
    }
}

void vtkCell::DeepCopy(vtkCell *c)
{
  this->Points->DeepCopy(c->Points);
  this->PointIds->DeepCopy(c->PointIds);
}

#define VTK_RIGHT 0
#define VTK_LEFT 1
#define VTK_MIDDLE 2

// Bounding box intersection modified from Graphics Gems Vol I. The method
// returns a non-zero value if the bounding box is hit. Origin[3] starts
// the ray, dir[3] is the vector components of the ray in the x-y-z
// directions, coord[3] is the location of hit, and t is the parametric
// coordinate along line. (Notes: the intersection ray dir[3] is NOT
// normalized.  Valid intersections will only occur between 0<=t<=1.)
char vtkCell::HitBBox (float bounds[6], float origin[3], float dir[3], 
                      float coord[3], float& t)
{
  char    inside=1;
  char    quadrant[3];
  int     i, whichPlane=0;
  float   maxT[3], candidatePlane[3];

  //  First find closest planes
  //
  for (i=0; i<3; i++) 
    {
    if ( origin[i] < bounds[2*i] ) 
      {
      quadrant[i] = VTK_LEFT;
      candidatePlane[i] = bounds[2*i];
      inside = 0;
      }
    else if ( origin[i] > bounds[2*i+1] ) 
      {
      quadrant[i] = VTK_RIGHT;
      candidatePlane[i] = bounds[2*i+1];
      inside = 0;
      }
    else 
      {
      quadrant[i] = VTK_MIDDLE;
      }
    }

  //  Check whether origin of ray is inside bbox
  //
  if (inside) 
    {
    coord[0] = origin[0];
    coord[1] = origin[1];
    coord[2] = origin[2];
    t = 0;
    return 1;
    }
  
  //  Calculate parametric distances to plane
  //
  for (i=0; i<3; i++)
    {
    if ( quadrant[i] != VTK_MIDDLE && dir[i] != 0.0 )
      {
      maxT[i] = (candidatePlane[i]-origin[i]) / dir[i];
      }
    else
      {
      maxT[i] = -1.0;
      }
    }

  //  Find the largest parametric value of intersection
  //
  for (i=0; i<3; i++)
    {
    if ( maxT[whichPlane] < maxT[i] )
      {
      whichPlane = i;
      }
    }

  //  Check for valid intersection along line
  //
  if ( maxT[whichPlane] > 1.0 || maxT[whichPlane] < 0.0 )
    {
    return 0;
    }
  else
    {
    t = maxT[whichPlane];
    }

  //  Intersection point along line is okay.  Check bbox.
  //
  for (i=0; i<3; i++) 
    {
    if (whichPlane != i) 
      {
      coord[i] = origin[i] + maxT[whichPlane]*dir[i];
      if ( coord[i] < bounds[2*i] || coord[i] > bounds[2*i+1] )
	{
        return 0;
	}
      } 
    else 
      {
      coord[i] = candidatePlane[i];
      }
    }

    return 1;
}

// Compute cell bounding box (xmin,xmax,ymin,ymax,zmin,zmax). Return pointer
// to array of six float values.
float *vtkCell::GetBounds ()
{
  float *x;
  int i, numPts=this->Points->GetNumberOfPoints();

  this->Bounds[0] = this->Bounds[2] = this->Bounds[4] =  VTK_LARGE_FLOAT;
  this->Bounds[1] = this->Bounds[3] = this->Bounds[5] = -VTK_LARGE_FLOAT;

  for (i=0; i<numPts; i++)
    {
    x = this->Points->GetPoint(i);

    this->Bounds[0] = (x[0] < this->Bounds[0] ? x[0] : this->Bounds[0]);
    this->Bounds[1] = (x[0] > this->Bounds[1] ? x[0] : this->Bounds[1]);
    this->Bounds[2] = (x[1] < this->Bounds[2] ? x[1] : this->Bounds[2]);
    this->Bounds[3] = (x[1] > this->Bounds[3] ? x[1] : this->Bounds[3]);
    this->Bounds[4] = (x[2] < this->Bounds[4] ? x[2] : this->Bounds[4]);
    this->Bounds[5] = (x[2] > this->Bounds[5] ? x[2] : this->Bounds[5]);
    
    }
  return this->Bounds;
}

// Compute cell bounding box (xmin,xmax,ymin,ymax,zmin,zmax). Copy result into
// user provided array.
void vtkCell::GetBounds(float bounds[6])
{
  this->GetBounds();
  for (int i=0; i < 6; i++)
    {
    bounds[i] = this->Bounds[i];
    }
}

// Compute Length squared of cell (i.e., bounding box diagonal squared).
float vtkCell::GetLength2 ()
{
  float diff, l=0.0;
  int i;

  this->GetBounds();
  for (i=0; i<3; i++)
    {
    diff = this->Bounds[2*i+1] - this->Bounds[2*i];
    l += diff * diff;
    }
 
  return l;
}

// Return center of the cell in parametric coordinates.
// Note that the parametric center is not always located 
// at (0.5,0.5,0.5). The return value is the subId that
// the center is in (if a composite cell). If you want the
// center in x-y-z space, invoke the EvaluateLocation() method.
int vtkCell::GetParametricCenter(float pcoords[3])
{
  pcoords[0] = pcoords[1] = pcoords[2] = 0.5;
  return 0;
}

void vtkCell::PrintSelf(ostream& os, vtkIndent indent)
{
  int numIds=this->PointIds->GetNumberOfIds();

  vtkObject::PrintSelf(os,indent);

  os << indent << "Number Of Points: " << numIds << "\n";

  if ( numIds > 0 )
    {
    float *bounds=this->GetBounds();

    os << indent << "Bounds: \n";
    os << indent << "  Xmin,Xmax: (" << bounds[0] << ", " << bounds[1] << ")\n";
    os << indent << "  Ymin,Ymax: (" << bounds[2] << ", " << bounds[3] << ")\n";
    os << indent << "  Zmin,Zmax: (" << bounds[4] << ", " << bounds[5] << ")\n";

    os << indent << "  Point ids are: ";
    for (int i=0; i < numIds; i++)
      {
      os << this->PointIds->GetId(i);
      if ( i && !(i % 12) )
	{
	os << "\n\t";
	}
      else
	{
	if ( i != (numIds-1) )
	  {
	  os << ", ";
	  }
	}
      }
    os << indent << "\n";
    }
}
