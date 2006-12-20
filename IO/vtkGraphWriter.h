/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkGraphWriter.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkGraphWriter - write vtk table data to a file
// .SECTION Description
// vtkGraphWriter is a source object that writes ASCII or binary 
// table data files in vtk format. See text for format details.
// .SECTION Caveats
// Binary files written on one system may not be readable on other systems.

#ifndef __vtkGraphWriter_h
#define __vtkGraphWriter_h

#include "vtkDataWriter.h"

class vtkGraph;

class VTK_IO_EXPORT vtkGraphWriter : public vtkDataWriter
{
public:
  static vtkGraphWriter *New();
  vtkTypeRevisionMacro(vtkGraphWriter,vtkDataWriter);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Get the input to this writer.
  vtkGraph* GetInput();
  vtkGraph* GetInput(int port);

protected:
  vtkGraphWriter() {};
  ~vtkGraphWriter() {};

  void WriteData();

  virtual int FillInputPortInformation(int port, vtkInformation *info);

private:
  vtkGraphWriter(const vtkGraphWriter&);  // Not implemented.
  void operator=(const vtkGraphWriter&);  // Not implemented.
};

#endif
