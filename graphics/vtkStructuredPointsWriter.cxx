/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkStructuredPointsWriter.cxx
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
#include "vtkStructuredPointsWriter.h"

// Specify the input data or filter.
void vtkStructuredPointsWriter::SetInput(vtkStructuredPoints *input)
{
  if ( this->Input != input )
    {
    vtkDebugMacro(<<" setting Input to " << (void *)input);
    if (this->Input) {this->Input->UnRegister(this);}
    this->Input = (vtkDataSet *) input;
    if (this->Input) {this->Input->Register(this);}
    this->Modified();
    }
}

void vtkStructuredPointsWriter::WriteData()
{
  FILE *fp;
  vtkStructuredPoints *input=(vtkStructuredPoints *)this->Input;
  int dim[3];
  float spacing[3], origin[3];

  vtkDebugMacro(<<"Writing vtk structured points...");

  if ( !(fp=this->OpenVTKFile()) || !this->WriteHeader(fp) )
      {
      return;
      }
//
// Write structured points specific stuff
//
  fprintf(fp,"DATASET STRUCTURED_POINTS\n");

  input->GetDimensions(dim);
  fprintf(fp,"DIMENSIONS %d %d %d\n", dim[0], dim[1], dim[2]);

  input->GetSpacing(spacing);
  fprintf(fp,"SPACING %g %g %g\n", spacing[0], spacing[1], spacing[2]);

  input->GetOrigin(origin);
  fprintf(fp,"ORIGIN %g %g %g\n", origin[0], origin[1], origin[2]);

  this->WriteCellData(fp, input);
  this->WritePointData(fp, input);

  this->CloseVTKFile(fp);
}

void vtkStructuredPointsWriter::PrintSelf(ostream& os, vtkIndent indent)
{
  vtkDataWriter::PrintSelf(os,indent);
}
