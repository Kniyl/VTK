/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkWriter.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkWriter.h"

#include "vtkCommand.h"
#include "vtkDataObject.h"
#include "vtkErrorCode.h"

vtkCxxRevisionMacro(vtkWriter, "1.40.10.1");

// Construct with no start and end write methods or arguments.
vtkWriter::vtkWriter()
{
}

vtkWriter::~vtkWriter()
{
}

vtkDataObject *vtkWriter::GetInput()
{
  if (this->NumberOfInputs < 1)
    {
    return NULL;
    }
  return this->Inputs[0];
}


// Write data to output. Method executes subclasses WriteData() method, as 
// well as StartMethod() and EndMethod() methods.
void vtkWriter::Write()
{
  this->SetErrorCode(vtkErrorCode::NoError);
  
  vtkDataObject *input = this->GetInput();
  int idx;

  // make sure input is available
  if ( !input )
    {
    vtkErrorMacro(<< "No input!");
    return;
    }

  for (idx = 0; idx < this->NumberOfInputs; ++idx)
    {
    if (this->Inputs[idx] != NULL)
      {
      this->Inputs[idx]->Update();
      }
    }

  unsigned long lastUpdateTime =  this->Inputs[0]->GetUpdateTime();
  for (idx = 1; idx < this->NumberOfInputs; ++idx)
    {
    unsigned long updateTime = this->Inputs[idx]->GetUpdateTime();
    if ( updateTime > lastUpdateTime )
      {
      lastUpdateTime = updateTime;
      }
    }

  if (lastUpdateTime < this->WriteTime && this->GetMTime() < this->WriteTime)
    {
    // we are up to date
    return;
    }

  this->InvokeEvent(vtkCommand::StartEvent,NULL);
  this->WriteData();
  this->InvokeEvent(vtkCommand::EndEvent,NULL);

  // Release any inputs if marked for release
  for (idx = 0; idx < this->NumberOfInputs; ++idx)
    {
    if (this->Inputs[idx] && this->Inputs[idx]->ShouldIReleaseData())
      {
      this->Inputs[idx]->ReleaseData();
      }
    }

  this->WriteTime.Modified();
}

// Convenient alias for Write() method.
void vtkWriter::Update()
{
  this->Write();
}

void vtkWriter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

}

void vtkWriter::EncodeArrayName(char* resname, const char* name)
{
  if ( !name || !resname )
    {   
    return;
    }
  int cc = 0;
  ostrstream str;

  char buffer[10];

  while( name[cc] )
    {
    // Encode spaces and %'s (and most non-printable ascii characters)
    // The reader does not support spaces in array names.
    if ( name[cc] < 33  || name[cc] > 126 ||
         name[cc] == '\"' || name[cc] == '%' )
      {
      sprintf(buffer, "%2X", name[cc]);
      str << "%%" << buffer; // Two % because it goes through printf format
      }
    else
      {
      str << name[cc];
      }
    cc++;
    }
  str << ends;
  strcpy(resname, str.str());
  str.rdbuf()->freeze(0);
}
