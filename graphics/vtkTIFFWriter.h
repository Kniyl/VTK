/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkTIFFWriter.h
  Language:  C++
  Date:      $Date$
  Version:   $Revision$


Copyright (c) 1993-1996 Ken Martin, Will Schroeder, Bill Lorensen.

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
// .NAME vtkTIFFWriter - write out structured points as a TIFF file
// .SECTION Description
// vtkTIFFWriter writes structured points as a non-compressed TIFF data file.
// The orientation of the image is with origin at lower left to correspond
// with vtk conventions. This can be changed with the SetOrientation method.

#ifndef __vtkTIFFWriter_h
#define __vtkTIFFWriter_h

#include <stdio.h>
#include "vtkWriter.h"
#include "vtkStructuredPoints.h"
#include "vtkImageToStructuredPoints.h"

class VTK_EXPORT vtkTIFFWriter : public vtkWriter
{
public:
  vtkTIFFWriter();
  ~vtkTIFFWriter();
  char *GetClassName() {return "vtkTIFFWriter";};
  void PrintSelf(ostream& os, vtkIndent indent);

  void SetInput(vtkStructuredPoints *input);
  void SetInput(vtkStructuredPoints &input) {this->SetInput(&input);};
  // Set color scalars ?
  void SetInput(vtkImageSource *cache)
    {this->SetInput(cache->GetImageToStructuredPoints()->GetOutput());}
  vtkStructuredPoints *GetInput() {return (vtkStructuredPoints *)this->Input;};

  // Description:
  // Specify name of file to write.
  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);
  void SetFilename(char *str){this->SetFileName(str);}
  char *GetFilename(){return this->GetFileName();}  

  // Description:
  // Specify Orientation of image. Default is 1.
  // 1	row 0 top, col 0 lhs
  // 2	row 0 top, col 0 rhs
  // 3	row 0 bottom, col 0 rhs
  // 4	row 0 bottom, col 0 lhs
  // 5	row 0 lhs, col 0 top
  // 6	row 0 rhs, col 0 top
  // 7	row 0 rhs, col 0 bottom
  // 8	row 0 lhs, col 0 bottom
  vtkSetClampMacro(Orientation,int,1,8);
  vtkGetMacro(Orientation,int);

protected:
  void WriteData();
  char *FileName;
  int Orientation;
};

#endif

