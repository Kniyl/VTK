/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkImageMagnify1D.h
  Language:  C++
  Date:      $Date$
  Version:   $Revision$
  Thanks:    Thanks to C. Charles Law who developed this class.

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
// .NAME vtkImageMagnify1D - Magnifies an image.
// .SECTION Description
// vtkImageMagnify1D maps each pixel of the input onto a n (integer)
// region of the output.  Location (0,0) remains in the same place.
// The filter can use pixel replication (nearest neighbor) or interpolation.


#ifndef __vtkImageMagnify1D_h
#define __vtkImageMagnify1D_h


#include "vtkImageFilter.h"

class VTK_EXPORT vtkImageMagnify1D : public vtkImageFilter
{
public:
  vtkImageMagnify1D();
  static vtkImageMagnify1D *New() {return new vtkImageMagnify1D;};
  const char *GetClassName() {return "vtkImageMagnify1D";};
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Set/Get the convolution axis.
  vtkSetMacro(MagnificationFactor,int);
  vtkGetMacro(MagnificationFactor,int);

  // Description:
  // Pixel replication or Interpolation.
  vtkSetMacro(Interpolate,int);
  vtkGetMacro(Interpolate,int);
  vtkBooleanMacro(Interpolate,int);
  
  // Description:
  // Specify which axis to magnify
  void SetFilteredAxis(int axis);
  int GetFilteredAxis() {return this->FilteredAxes[0];}
  
protected:
  int MagnificationFactor;
  int Interpolate;

  void ExecuteImageInformation(vtkImageCache *in, vtkImageCache *out);
  void ComputeRequiredInputUpdateExtent(vtkImageCache *out, vtkImageCache *in);
  void Execute(vtkImageRegion *inRegion, vtkImageRegion *outRegion);  
};

#endif



