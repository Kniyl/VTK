/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkColorScalars.cc
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
#include "vtkColorScalars.hh"
#include "vtkLookupTable.hh"

// Description:
// Convert internal color representation into scalar value. Uses luminance 
// equation (see text).
float vtkColorScalars::GetScalar(int i)
{
  unsigned char *rgba;
  float s;

  rgba  = this->GetColor(i);
  s = (rgba[3]/255.0) * (0.30*rgba[0] + 0.59*rgba[1] + 0.11*rgba[2]);
  if ( s > 1.0 ) s = 1.0;
  return s;
}

// Description:
// Map through lookup table to set the color.
void vtkColorScalars::SetScalar(int i, float s)
{
  if ( this->LookupTable == NULL ) this->CreateDefaultLookupTable();

  this->SetColor(i,this->LookupTable->MapValue(s));
}

// Description:
// Map through lookup table  to set the color.
void vtkColorScalars::InsertScalar(int i, float s)
{
  if ( this->LookupTable == NULL ) this->CreateDefaultLookupTable();

  this->InsertColor(i,this->LookupTable->MapValue(s));
}

// Description:
// Map through lookup table to set the color.
int vtkColorScalars::InsertNextScalar(float s)
{
  if ( this->LookupTable == NULL ) this->CreateDefaultLookupTable();

  return this->InsertNextColor(this->LookupTable->MapValue(s));
}

// Description:
// Given list of point id's, return colors for each point.
void vtkColorScalars::GetColors(vtkIdList& ptId, vtkAPixmap& p)
{
  for (int i=0; i<ptId.GetNumberOfIds(); i++)
    {
    p.InsertColor(i,this->GetColor(ptId.GetId(i)));
    }
}

// Description:
// Compute range of color rgba data (rmin,rmax, gmin,gmax, bmin,bmax, 
// amin,amax). Return pointer to array of length 8.
unsigned char *vtkColorScalars::GetComponentRange ()
{
  unsigned char *rgba;
  int i, j;
  static unsigned char range[8];

  range[0] = range[2] = range[4] = range[6] = 255;
  range[1] = range[3] = range[5] = range[7] = 0;

  for (i=0; i<this->GetNumberOfColors(); i++)
    {
    rgba = this->GetColor(i);
    for (j=0; j<4; j++)
      {
      if ( rgba[j] < range[2*j] ) range[2*j] = rgba[j];
      if ( rgba[j] > range[2*j+1] ) range[2*j+1] = rgba[j];
      }
    }
  return range;
}

// Description:
// Compute range of color rgba data (rmin,rmax, gmin,gmax, bmin,bmax, 
// amin,amax). Copy result into user provided array.
void vtkColorScalars::GetComponentRange(unsigned char range[8])
{
  unsigned char *r=this->GetComponentRange();
  for (int i=0; i < 8; i++) range[i] = r[i];
}

