/*=========================================================================

Program:   Visualization Toolkit
Module:    vtkHierarchicalDataSetAlgorithm.cxx

Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
All rights reserved.
See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkHierarchicalDataSetAlgorithm.h"

#include "vtkCommand.h"
#include "vtkCompositeDataPipeline.h"
#include "vtkDataSet.h"
#include "vtkHierarchicalDataSet.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"

vtkCxxRevisionMacro(vtkHierarchicalDataSetAlgorithm, "1.7");
vtkStandardNewMacro(vtkHierarchicalDataSetAlgorithm);

//----------------------------------------------------------------------------
// Instantiate object so that cell data is not passed to output.
vtkHierarchicalDataSetAlgorithm::vtkHierarchicalDataSetAlgorithm()
{
}

//----------------------------------------------------------------------------
vtkHierarchicalDataSet* vtkHierarchicalDataSetAlgorithm::GetOutput()
{
  return this->GetOutput(0);
}

//----------------------------------------------------------------------------
vtkHierarchicalDataSet* vtkHierarchicalDataSetAlgorithm::GetOutput(int port)
{
  vtkDataObject* output = 
    vtkCompositeDataPipeline::SafeDownCast(this->GetExecutive())->
    GetCompositeOutputData(port);
  return vtkHierarchicalDataSet::SafeDownCast(output);
}

//----------------------------------------------------------------------------
int vtkHierarchicalDataSetAlgorithm::FillOutputPortInformation(
  int vtkNotUsed(port), vtkInformation* info)
{
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataObject");
  info->Set(vtkCompositeDataPipeline::COMPOSITE_DATA_TYPE_NAME(), 
            "vtkHierarchicalDataSet");
  return 1;
}

//----------------------------------------------------------------------------
int vtkHierarchicalDataSetAlgorithm::FillInputPortInformation(
  int vtkNotUsed(port), vtkInformation* info)
{
  // now add our info
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataObject");
  info->Set(vtkCompositeDataPipeline::INPUT_REQUIRED_COMPOSITE_DATA_TYPE(), 
            "vtkHierarchicalDataSet");
  return 1;
}

//----------------------------------------------------------------------------
void vtkHierarchicalDataSetAlgorithm::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}
