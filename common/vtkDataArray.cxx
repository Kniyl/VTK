/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkDataArray.cxx
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
#include "vtkDataArray.h"
#include "vtkFloatArray.h"

// Construct object with default tuple dimension (number of components) of 1.
vtkDataArray::vtkDataArray(int numComp)
{
  this->Size = 0;
  this->MaxId = -1;
  this->Extend = 1000;

  this->NumberOfComponents = (numComp < 1 ? 1 : numComp);
}

void vtkDataArray::DeepCopy(vtkDataArray *da)
{
  if ( this != da )
    {
    int numTuples = da->GetNumberOfTuples();
    this->NumberOfComponents = da->NumberOfComponents;
    this->SetNumberOfTuples(numTuples);

    for (int i=0; i < numTuples; i++)
      {
      this->SetTuple(i, da->GetTuple(i));
      }
    }
}

// These can be overridden for more efficiency
float vtkDataArray::GetComponent(const int i, const int j)
{
  float *tuple=new float[this->NumberOfComponents], c;

  this->GetTuple(i,tuple);
  c =  tuple[j];
  delete [] tuple;

  return c;
}

void vtkDataArray::SetComponent(const int i, const int j, const float c)
{
  float *tuple=new float[this->NumberOfComponents];

  if ( i < this->GetNumberOfTuples() )
    {
    this->GetTuple(i,tuple);
    }
  else
    {
    for (int k=0; k<this->NumberOfComponents; k++)
      {
      tuple[k] = 0.0;
      }
    }

  tuple[j] = c;
  this->SetTuple(i,tuple);

  delete [] tuple;
}

void vtkDataArray::InsertComponent(const int i, const int j, const float c)
{
  float *tuple=new float[this->NumberOfComponents];

  if ( i < this->GetNumberOfTuples() )
    {
    this->GetTuple(i,tuple);
    }
  else
    {
    for (int k=0; k<this->NumberOfComponents; k++)
      {
      tuple[k] = 0.0;
      }
    }

  tuple[j] = c;
  this->InsertTuple(i,tuple);

  delete [] tuple;
}

void vtkDataArray::GetData(int tupleMin, int tupleMax, int compMin, int compMax, 
			   vtkFloatArray &data)
{
  int i, j;
  int numComp=this->GetNumberOfComponents();
  float *tuple=new float[numComp];
  float *ptr=data.WritePointer(0,(tupleMax-tupleMin+1)*(compMax-compMin+1));
  
  for (j=tupleMin; j <= tupleMax; j++)
    {
    this->GetTuple(j,tuple);
    for (i=compMin; i <= compMax; i++)
      {
      *ptr++ = tuple[i];
      }
    }
  delete [] tuple;
}

// default double behaviour
void vtkDataArray::GetTuple(const int i, double * tuple)
{
  int c;
  int numComp=this->GetNumberOfComponents();
  float *ftuple=new float[numComp];
  this->GetTuple(i,ftuple);
  for (c = 0; c < numComp;  c++)
    {
    tuple[c] = ftuple[c];
    }
  delete [] ftuple;
}

void vtkDataArray::SetTuple(const int i, const double * tuple)
{
  int c;
  int numComp=this->GetNumberOfComponents();
  float *ftuple=new float[numComp];
  for (c = 0; c < numComp;  c++)
    {
    ftuple[c] = (float)(tuple[c]);
    }
  this->SetTuple(i,ftuple);
  delete [] ftuple;
}

void vtkDataArray::InsertTuple(const int i, const double * tuple)
{
  int c;
  int numComp=this->GetNumberOfComponents();
  float *ftuple=new float[numComp];
  for (c = 0; c < numComp;  c++)
    {
    ftuple[c] = (float)(tuple[c]);
    }
  this->InsertTuple(i,ftuple);
  delete [] ftuple;
}

int vtkDataArray::InsertNextTuple(const double * tuple)
{
  int c;
  int numComp=this->GetNumberOfComponents();
  float *ftuple=new float[numComp];
  for (c = 0; c < numComp;  c++)
    {
    ftuple[c] = (float)(tuple[c]);
    }
  int ret = this->InsertNextTuple(ftuple);
  delete [] ftuple;
  return ret;
}

unsigned long vtkDataArray::GetActualMemorySize()
{
  unsigned long numPrims;
  float size;
  numPrims = this->GetNumberOfTuples() * this->GetNumberOfComponents();

  switch (this->GetDataType())
    {
    case VTK_BIT:
      size = (float)sizeof(char)/8.0;
      break;

    case VTK_CHAR:
      size = (float)sizeof(char);
      break;

    case VTK_UNSIGNED_CHAR:
      size = (float)sizeof(unsigned char);
      break;

    case VTK_SHORT:
      size = (float)sizeof(short);
      break;

    case VTK_UNSIGNED_SHORT:
      size = (float)sizeof(unsigned short);
      break;

    case VTK_INT:
      size = (float)sizeof(int);
      break;

    case VTK_UNSIGNED_INT:
      size = (float)sizeof(unsigned int);
      break;

    case VTK_LONG:
      size = (float)sizeof(long);
      break;

    case VTK_UNSIGNED_LONG:
      size = (float)sizeof(unsigned long);
      break;

    case VTK_FLOAT:
      size = (float)sizeof(float);
      break;

    case VTK_DOUBLE:
      size = (float)sizeof(double);
      break;

    default:
      vtkErrorMacro(<<"Unsupported data type!");
    }

  return (unsigned long)ceil((size * numPrims)/1000.0); //kilobytes
}

void vtkDataArray::PrintSelf(ostream& os, vtkIndent indent)
{
  vtkObject::PrintSelf(os,indent);

  os << indent << "Number Of Components: " << this->NumberOfComponents << "\n";
  os << indent << "Number Of Tuples: " << this->GetNumberOfTuples() << "\n";
  os << indent << "Size: " << this->Size << "\n";
  os << indent << "MaxId: " << this->MaxId << "\n";
  os << indent << "Extend size: " << this->Extend << "\n";
}
