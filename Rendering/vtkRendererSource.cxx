/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkRendererSource.cxx
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

  Copyright (c) 1993-2002 Ken Martin, Will Schroeder, Bill Lorensen 
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkRendererSource.h"

#include "vtkFloatArray.h"
#include "vtkImageData.h"
#include "vtkMapper.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkRenderWindow.h"
#include "vtkRenderer.h"
#include "vtkUnsignedCharArray.h"

vtkCxxRevisionMacro(vtkRendererSource, "1.49");
vtkStandardNewMacro(vtkRendererSource);

vtkCxxSetObjectMacro(vtkRendererSource,Input,vtkRenderer);

vtkRendererSource::vtkRendererSource()
{
  this->Input = NULL;
  this->WholeWindow = 0;
  this->RenderFlag = 0;
  this->DepthValues = 0;
  this->DepthValuesInScalars = 0;
}


vtkRendererSource::~vtkRendererSource()
{
  if (this->Input)
    {
    this->Input->UnRegister(this);
    this->Input = NULL;
    }
}

void vtkRendererSource::ExecuteData(vtkDataObject *outp)
{
  int numOutPts;
  float x1,y1,x2,y2;
  unsigned char *pixels, *ptr;
  int dims[3];
  vtkImageData *output = this->AllocateOutputData(outp);
  vtkUnsignedCharArray *outScalars = 
    vtkUnsignedCharArray::SafeDownCast(output->GetPointData()->GetScalars());
  if (this->DepthValuesInScalars)
    {
    outScalars->SetName("RGBValues");
    }
  else
    {
    outScalars->SetName("RGBZValues");
    }
  vtkRenderWindow *renWin;
  
  vtkDebugMacro(<<"Converting points");

  if (this->Input == NULL )
    {
    vtkErrorMacro(<<"Please specify a renderer as input!");
    return;
    }

  renWin = this->Input->GetRenderWindow();
  if (renWin == NULL)
    {
    return;
    }
  
  if (this->RenderFlag)
    {
    renWin->Render();
    }
  
  // calc the pixel range for the renderer
  x1 = this->Input->GetViewport()[0]*
    ((this->Input->GetRenderWindow())->GetSize()[0] - 1);
  y1 = this->Input->GetViewport()[1]*
    ((this->Input->GetRenderWindow())->GetSize()[1] - 1);
  x2 = this->Input->GetViewport()[2]*
    ((this->Input->GetRenderWindow())->GetSize()[0] - 1);
  y2 = this->Input->GetViewport()[3]*
    ((this->Input->GetRenderWindow())->GetSize()[1] - 1);

  if (this->WholeWindow)
    {
    x1 = 0;
    y1 = 0;
    x2 = (this->Input->GetRenderWindow())->GetSize()[0] - 1;
    y2 = (this->Input->GetRenderWindow())->GetSize()[1] - 1;
    }
  
  // Get origin, aspect ratio and dimensions from this->Input
  dims[0] = (int)(x2 - x1 + 1); dims[1] = (int)(y2 -y1 + 1); dims[2] = 1;
  output->SetDimensions(dims);
  output->SetSpacing(1,1,1);
  output->SetOrigin(0,0,0);

  // Allocate data.  Scalar type is FloatScalars.
  numOutPts = dims[0] * dims[1];

  pixels = (this->Input->GetRenderWindow())->GetPixelData((int)x1,(int)y1,
                                                          (int)x2,(int)y2,1);

  // allocate scalars
  int nb_comp = output->GetNumberOfScalarComponents();
  ptr = outScalars->WritePointer(0, numOutPts * nb_comp);

  // copy scalars over (if only RGB is requested, use the pixels directly)
  if (!this->DepthValuesInScalars)
    {
    memcpy(ptr, pixels, numOutPts * nb_comp);
    }
  
  // Lets get the ZBuffer also, if requested.
  if (this->DepthValues || this->DepthValuesInScalars)
    {
    float *zBuf = (this->Input->GetRenderWindow())->GetZbufferData(
      (int)x1,(int)y1, (int)x2,(int)y2);

    // If RGBZ is requested, intermix RGB with shift/scaled Z
    if (this->DepthValuesInScalars)
      {
      float *zptr = zBuf, *zptr_end = zBuf + numOutPts;
      float min = *zBuf, max = *zBuf;
      while (zptr < zptr_end)
        {
        if (min < *zptr) { min = *zptr; }
        if (max > *zptr) { max = *zptr; }
        zptr++;
        }
      float scale = 255.0 / (max - min);

      zptr = zBuf;
      unsigned char *ppixels = pixels;
      while (zptr < zptr_end)
        {
        *ptr++ = *ppixels++;
        *ptr++ = *ppixels++;
        *ptr++ = *ppixels++;
        *ptr++ = static_cast<unsigned char>((*zptr++ - min) * scale);
        }
      }

    // If Z is requested as independent array, create it
    if (this->DepthValues)
      {
      vtkFloatArray *zArray = vtkFloatArray::New();
      zArray->Allocate(numOutPts);
      zArray->SetNumberOfTuples(numOutPts);
      float *zPtr = zArray->WritePointer(0, numOutPts);
      memcpy(zPtr,zBuf,numOutPts*sizeof(float));
      zArray->SetName("ZBuffer");
      output->GetPointData()->AddArray(zArray);
      zArray->Delete();
      }

    delete [] zBuf;
    }

  delete [] pixels;
}

void vtkRendererSource::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "RenderFlag: " << (this->RenderFlag ? "On\n" : "Off\n");

  if ( this->Input )
    {
    os << indent << "Input:\n";
    this->Input->PrintSelf(os,indent.GetNextIndent());
    }
  else
    {
    os << indent << "Input: (none)\n";
    }

  os << indent << "Whole Window: " << (this->WholeWindow ? "On\n" : "Off\n");
  os << indent << "Depth Values: " << (this->DepthValues ? "On\n" : "Off\n");
  os << indent << "Depth Values In Scalars: " << (this->DepthValuesInScalars ? "On\n" : "Off\n");
}


unsigned long vtkRendererSource::GetMTime()
{
  unsigned long mTime=this->MTime.GetMTime();
  unsigned long transMTime;

  if ( this->Input )
    {
    transMTime = this->Input->GetMTime();
    mTime = ( transMTime > mTime ? transMTime : mTime );
    }

  return mTime;
}

//----------------------------------------------------------------------------
// Consider renderer for PiplineMTime
void vtkRendererSource::UpdateInformation()
{
  unsigned long t1, t2;
  vtkImageData *output = this->GetOutput();
  vtkRenderer *ren = this->GetInput();
  vtkActorCollection *actors;
  vtkActor *actor;
  vtkMapper *mapper;
  vtkDataObject *data=0;
  float x1,y1,x2,y2;
  
  if (output == NULL || ren == NULL || ren->GetRenderWindow() == NULL)
    {
    return;
    }
  
  // calc the pixel range for the renderer
  x1 = ren->GetViewport()[0] * ((ren->GetRenderWindow())->GetSize()[0] - 1);
  y1 = ren->GetViewport()[1] * ((ren->GetRenderWindow())->GetSize()[1] - 1);
  x2 = ren->GetViewport()[2] * ((ren->GetRenderWindow())->GetSize()[0] - 1);
  y2 = ren->GetViewport()[3] *((ren->GetRenderWindow())->GetSize()[1] - 1);
  if (this->WholeWindow)
    {
    x1 = 0;
    y1 = 0;
    x2 = (this->Input->GetRenderWindow())->GetSize()[0] - 1;
    y2 = (this->Input->GetRenderWindow())->GetSize()[1] - 1;
    }    
  output->SetWholeExtent(0, static_cast<int>(x2-x1), 
                         0, static_cast<int>(y2-y1), 0, 0);
  output->SetScalarType(VTK_UNSIGNED_CHAR);
  output->SetNumberOfScalarComponents(3 + (this->DepthValuesInScalars ? 1:0));
  
  // Update information on the input and
  // compute information that is general to vtkDataObject.
  t1 = this->GetMTime();
  t2 = ren->GetMTime();
  if (t2 > t1)
    {
    t1 = t2;
    }
  actors = ren->GetActors();
  actors->InitTraversal();
  while ( (actor = actors->GetNextItem()) )
    {
    t2 = actor->GetMTime();
    if (t2 > t1)
      {
      t1 = t2;
      }
    mapper = actor->GetMapper();
    if (mapper)
      {
      t2 = mapper->GetMTime();
      if (t2 > t1)
        {
        t1 = t2;
        }
      data = mapper->GetInput();
      if (data)
        {
        data->UpdateInformation();
        }
      t2 = data->GetMTime();
      if (t2 > t1)
        {
        t1 = t2;
        }
      t2 = data->GetPipelineMTime();
      if (t2 > t1)
        {
        t1 = t2;
        }
      }  
    }

  output->SetPipelineMTime(t1);
  this->InformationTime.Modified();
}


