/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkExtractGrid.cxx
  Language:  C++
  Date:      $Date$
  Version:   $Revision$


Copyright (c) 1993-2000 Ken Martin, Will Schroeder, Bill Lorensen 
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.

 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

 * Neither name of Ken Martin, Will Schroeder, or Bill Lorensen nor the names
   of any contributors may be used to endorse or promote products derived
   from this software without specific prior written permission.

 * Modified source versions must be plainly marked as such, and must not be
   misrepresented as being the original software.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS ``AS IS''
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=========================================================================*/
#include "vtkExtractGrid.h"
#include "vtkObjectFactory.h"



//------------------------------------------------------------------------------
vtkExtractGrid* vtkExtractGrid::New()
{
  // First try to create the object from the vtkObjectFactory
  vtkObject* ret = vtkObjectFactory::CreateInstance("vtkExtractGrid");
  if(ret)
    {
    return (vtkExtractGrid*)ret;
    }
  // If the factory was unable to create the object, then create it here.
  return new vtkExtractGrid;
}




// Construct object to extract all of the input data.
vtkExtractGrid::vtkExtractGrid()
{
  this->VOI[0] = this->VOI[2] = this->VOI[4] = 0;
  this->VOI[1] = this->VOI[3] = this->VOI[5] = VTK_LARGE_INTEGER;

  this->SampleRate[0] = this->SampleRate[1] = this->SampleRate[2] = 1;
}

void vtkExtractGrid::Execute()
{
  vtkStructuredGrid *input= this->GetInput();
  vtkPointData *pd=input->GetPointData();
  vtkCellData *cd=input->GetCellData();
  vtkStructuredGrid *output= this->GetOutput();
  vtkPointData *outPD=output->GetPointData();
  vtkCellData *outCD=output->GetCellData();
  int i, j, k, dims[3], outDims[3], voi[6], dim, idx, newIdx;
  int newCellId;
  int sliceSize, outSize, jOffset, kOffset, rate[3];
  vtkPoints *newPts, *inPts;

  vtkDebugMacro(<< "Extracting Grid");
//
// Check Grid and clamp as necessary. Compute output parameters,
//
  inPts = input->GetPoints();
  input->GetDimensions(dims);

  for ( i=0; i < 6; i++ )
    {
    voi[i] = this->VOI[i];
    }

  for ( outSize=1, dim=0, i=0; i < 3; i++ )
    {
    if ( voi[2*i+1] >= dims[i] )
      {
      voi[2*i+1] = dims[i] - 1;
      }
    else if ( voi[2*i+1] < 0 )
      {
      voi[2*i+1] = 0;
      }

    if ( voi[2*i] > voi[2*i+1] )
      {
      voi[2*i] = voi[2*i+1];
      }
    else if ( voi[2*i] < 0 )
      {
      voi[2*i] = 0;
      }

    if ( (voi[2*i+1]-voi[2*i]) > 0 )
      {
      dim++;
      }

    if ( (rate[i] = this->SampleRate[i]) < 1 )
      {
      rate[i] = 1;
      }

    outDims[i] = (voi[2*i+1] - voi[2*i]) / rate[i] + 1;
    if ( outDims[i] < 1 )
      {
      outDims[i] = 1;
      }

    outSize *= outDims[i];
    }

  output->SetDimensions(outDims);
//
// If output same as input, just pass data through
//
  if ( outDims[0] == dims[0] && outDims[1] == dims[1] && outDims[2] == dims[2] &&
  rate[0] == 1 && rate[1] == 1 && rate[2] == 1 )
    {
    output->SetPoints(inPts);
    output->GetPointData()->PassData(input->GetPointData());
    output->GetCellData()->PassData(input->GetCellData());
    vtkDebugMacro(<<"Passed data through bacause input and output are the same");
    return;
    }
//
// Allocate necessary objects
//
  newPts = (vtkPoints *) inPts->MakeObject(); newPts->SetNumberOfPoints(outSize);
  outPD->CopyAllocate(pd,outSize,outSize);
  outCD->CopyAllocate(cd,outSize,outSize);
//
// Traverse input data and copy point attributes to output
//
  sliceSize = dims[0]*dims[1];
  newIdx = 0;
  for ( k=voi[4]; k <= voi[5]; k += rate[2] )
    {
    kOffset = k * sliceSize;
    for ( j=voi[2]; j <= voi[3]; j += rate[1] )
      {
      jOffset = j * dims[0];
      for ( i=voi[0]; i <= voi[1]; i += rate[0] )
        {
        idx = i + jOffset + kOffset;
        newPts->SetPoint(newIdx,inPts->GetPoint(idx));
        outPD->CopyData(pd, idx, newIdx++);
        }
      }
    }

//
// Traverse input data and copy cell attributes to output
//
  newCellId = 0;
  sliceSize = (dims[0]-1)*(dims[1]-1);
  for ( k=voi[4]; k < voi[5]; k += rate[2] )
    {
    kOffset = k * sliceSize;
    for ( j=voi[2]; j < voi[3]; j += rate[1] )
      {
      jOffset = j * (dims[0] - 1);
      for ( i=voi[0]; i < voi[1]; i += rate[0] )
        {
        idx = i + jOffset + kOffset;
        outCD->CopyData(cd, idx, newCellId++);
        }
      }
    }

  vtkDebugMacro(<<"Extracted " << newIdx << " point attributes on "
                << dim << "-D dataset\n\tDimensions are (" << outDims[0]
                << "," << outDims[1] << "," << outDims[2] <<")");

  output->SetPoints(newPts);
  newPts->Delete();
}


void vtkExtractGrid::PrintSelf(ostream& os, vtkIndent indent)
{
  vtkStructuredGridToStructuredGridFilter::PrintSelf(os,indent);

  os << indent << "VOI: \n";
  os << indent << "  Imin,Imax: (" << this->VOI[0] << ", " << this->VOI[1] << ")\n";
  os << indent << "  Jmin,Jmax: (" << this->VOI[2] << ", " << this->VOI[3] << ")\n";
  os << indent << "  Kmin,Kmax: (" << this->VOI[4] << ", " << this->VOI[5] << ")\n";

  os << indent << "Sample Rate: (" << this->SampleRate[0] << ", "
               << this->SampleRate[1] << ", "
               << this->SampleRate[2] << ")\n";
}


