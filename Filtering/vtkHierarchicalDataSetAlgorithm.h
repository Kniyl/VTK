/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkHierarchicalDataSetAlgorithm.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkHierarchicalDataSetAlgorithm - Superclass for algorithms that produce only vtkHierarchicalDataSet as output
// .SECTION Description
// Algorithms that take any type of data object (including composite dataset)
// and produce a vtkHierarchicalDataSet in the output can subclass from this
// class.


#ifndef __vtkHierarchicalDataSetAlgorithm_h
#define __vtkHierarchicalDataSetAlgorithm_h

#include "vtkMultiGroupDataSetAlgorithm.h"

class vtkHierarchicalDataSet;

class VTK_FILTERING_EXPORT vtkHierarchicalDataSetAlgorithm : public vtkMultiGroupDataSetAlgorithm
{
public:
  static vtkHierarchicalDataSetAlgorithm *New();
  vtkTypeRevisionMacro(vtkHierarchicalDataSetAlgorithm,vtkMultiGroupDataSetAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Get the output data object for a port on this algorithm.
  vtkHierarchicalDataSet* GetOutput();
  vtkHierarchicalDataSet* GetOutput(int);

protected:
  vtkHierarchicalDataSetAlgorithm();
  ~vtkHierarchicalDataSetAlgorithm() {};

  // see algorithm for more info
  virtual int FillOutputPortInformation(int port, vtkInformation* info);
  virtual int FillInputPortInformation(int port, vtkInformation* info);

private:
  vtkHierarchicalDataSetAlgorithm(const vtkHierarchicalDataSetAlgorithm&);  // Not implemented.
  void operator=(const vtkHierarchicalDataSetAlgorithm&);  // Not implemented.
};

#endif


