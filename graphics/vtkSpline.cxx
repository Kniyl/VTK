/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkSpline.cxx
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

#include "vtkSpline.h"

// Construct a spline wth the folloing defaults:
// ClampValueOff
vtkSpline::vtkSpline ()
{
  this->ComputeTime = 0;
  this->ClampValue = 0;
  this->PiecewiseFunction = vtkPiecewiseFunction::New();
  this->Intervals = NULL;
  this->Coefficients = NULL;
  this->LeftConstraint = 1;
  this->LeftValue = 0.0;
  this->RightConstraint = 1;
  this->RightValue = 0.0;
  this->Closed = 0;
}

vtkSpline::~vtkSpline ()
{
  if (this->PiecewiseFunction)
    {
    this->PiecewiseFunction->Delete();
    }
  if (this->Coefficients)
    {
    delete [] this->Coefficients;
    }
  if (this->Intervals)
    {
    delete [] this->Intervals;
    }
}

// Add a point to the Piecewise Functions containing the data
void vtkSpline::AddPoint (float t, float x)
{
  this->PiecewiseFunction->AddPoint (t, x);
}

// Remove a point from the Piecewise Functions.
void vtkSpline::RemovePoint (float t)
{
  this->PiecewiseFunction->RemovePoint (t);
}

// Remove all points from the Piecewise Functions.
void vtkSpline::RemoveAllPoints ()
{
  this->PiecewiseFunction->RemoveAllPoints ();
}

// Overload standard modified time function. If data is modified,
// then this object is modified as well.
unsigned long vtkSpline::GetMTime()
{
  unsigned long mTime=this->vtkObject::GetMTime();
  unsigned long DataMTime;

  if ( this->PiecewiseFunction != NULL )
    {
    DataMTime = this->PiecewiseFunction->GetMTime();
    mTime = ( DataMTime > mTime ? DataMTime : mTime );
    }

  return mTime;
}

void vtkSpline::PrintSelf(ostream& os, vtkIndent indent)
{
  vtkObject::PrintSelf(os,indent);

  os << indent << "Clamp Value: " << (this->ClampValue ? "On\n" : "Off\n");
  os << indent << "Left Constraint: " << this->LeftConstraint << "\n";
  os << indent << "Right Constraint: " << this->RightConstraint << "\n";
  os << indent << "Left Value: " << this->LeftValue << "\n";
  os << indent << "Right Value: " << this->RightValue << "\n";
  os << indent << "Closed: " << (this->Closed ? "On\n" : "Off\n");

  os << indent << "Piecewise Function:\n";
  this->PiecewiseFunction->PrintSelf(os,indent.GetNextIndent());

  os << indent << "Closed: " << (this->Closed ? "On\n" : "Off\n");
}
