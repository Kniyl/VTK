/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkXMLImageDataReader.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkXMLImageDataReader.h"

#include "vtkDataArray.h"
#include "vtkImageData.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkXMLDataElement.h"
#include "vtkInformation.h"
#include "vtkStreamingDemandDrivenPipeline.h"

vtkCxxRevisionMacro(vtkXMLImageDataReader, "1.5");
vtkStandardNewMacro(vtkXMLImageDataReader);

//----------------------------------------------------------------------------
vtkXMLImageDataReader::vtkXMLImageDataReader()
{
  vtkImageData *output = vtkImageData::New();
  this->SetOutput(output);
  // Releasing data for pipeline parallism.
  // Filters will know it is empty. 
  output->ReleaseData();
  output->Delete();
}

//----------------------------------------------------------------------------
vtkXMLImageDataReader::~vtkXMLImageDataReader()
{
}

//----------------------------------------------------------------------------
void vtkXMLImageDataReader::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

//----------------------------------------------------------------------------
void vtkXMLImageDataReader::SetOutput(vtkImageData *output)
{
  this->GetExecutive()->SetOutputData(0, output);
}

//----------------------------------------------------------------------------
vtkImageData* vtkXMLImageDataReader::GetOutput()
{
  return this->GetOutput(0);
}

//----------------------------------------------------------------------------
vtkImageData* vtkXMLImageDataReader::GetOutput(int idx)
{
  return vtkImageData::SafeDownCast( this->GetOutputDataObject(idx) );
}


//----------------------------------------------------------------------------
const char* vtkXMLImageDataReader::GetDataSetName()
{
  return "ImageData";
}

//----------------------------------------------------------------------------
void vtkXMLImageDataReader::SetOutputExtent(int* extent)
{
  this->GetOutput()->SetExtent(extent);
}

//----------------------------------------------------------------------------
int vtkXMLImageDataReader::ReadPrimaryElement(vtkXMLDataElement* ePrimary)
{
  if(!this->Superclass::ReadPrimaryElement(ePrimary)) { return 0; }
  
  // Get the image's origin.
  if(ePrimary->GetVectorAttribute("Origin", 3, this->Origin) != 3)
    {
    this->Origin[0] = 0;
    this->Origin[1] = 0;
    this->Origin[2] = 0;
    }
  
  // Get the image's spacing.
  if(ePrimary->GetVectorAttribute("Spacing", 3, this->Spacing) != 3)
    {
    this->Spacing[0] = 1;
    this->Spacing[1] = 1;
    this->Spacing[2] = 1;
    }
  
  return 1;
}

//----------------------------------------------------------------------------
void vtkXMLImageDataReader::SetupOutputInformation(vtkInformation *outInfo)
{
  this->Superclass::SetupOutputInformation(outInfo);

  outInfo->Set(vtkDataObject::ORIGIN(), this->Origin, 3);
  outInfo->Set(vtkDataObject::SPACING(), this->Spacing, 3);
    
  // Backward-compatability support for scalar information in output.
  if (this->PointDataElements[0])
    {
    int components, dataType, i;
    for(i = 0; i < this->PointDataElements[0]->GetNumberOfNestedElements(); i++)
      {
      vtkXMLDataElement* eNested = this->PointDataElements[0]->GetNestedElement(i);
      if ( eNested->GetAttribute( "Scalars" ) )
        {
        if(!eNested->GetWordTypeAttribute("type", dataType))
          {
          this->InformationError = 1;
          return;
          }
        if(!eNested->GetScalarAttribute("NumberOfComponents", components))
          {
          this->InformationError = 1;
          return;
          }
        outInfo->Set(vtkDataObject::SCALAR_TYPE(), dataType);
        outInfo->Set(vtkDataObject::SCALAR_NUMBER_OF_COMPONENTS(), components);
        break;
        }
      }
    }
}


//----------------------------------------------------------------------------
int vtkXMLImageDataReader::FillOutputPortInformation(int, vtkInformation* info)
{
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkImageData");
  return 1;
}

