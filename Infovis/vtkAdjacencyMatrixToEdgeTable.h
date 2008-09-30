/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkAdjacencyMatrixToEdgeTable.h
  
-------------------------------------------------------------------------
  Copyright 2008 Sandia Corporation.
  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
  the U.S. Government retains certain rights in this software.
-------------------------------------------------------------------------

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __vtkAdjacencyMatrixToEdgeTable_h
#define __vtkAdjacencyMatrixToEdgeTable_h

#include "vtkTableAlgorithm.h"

// .NAME vtkAdjacencyMatrixToEdgeTable

// .SECTION Description
// Treats a dense 2-way array of doubles as an adacency matrix and converts it into a
// vtkTable suitable for use as an edge table with vtkTableToGraph.

// .SECTION Thanks
// Developed by Timothy M. Shead (tshead@sandia.gov) at Sandia National Laboratories.

class VTK_INFOVIS_EXPORT vtkAdjacencyMatrixToEdgeTable : public vtkTableAlgorithm
{
public:
  static vtkAdjacencyMatrixToEdgeTable* New();
  vtkTypeRevisionMacro(vtkAdjacencyMatrixToEdgeTable, vtkTableAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Specifies whether rows or columns become the "source" in the output edge table.
  // 0 = rows, 1 = columns.
  vtkGetMacro(SourceDimension, vtkIdType);
  vtkSetMacro(SourceDimension, vtkIdType);

  // Description:
  // Controls the name of the output table column that contains edge weights.
  vtkGetStringMacro(ValueArrayName);
  vtkSetStringMacro(ValueArrayName);

  // Description:
  // Specifies the minimum number of adjacent edges to include for each source vertex.
  vtkGetMacro(MinimumCount, vtkIdType);
  vtkSetMacro(MinimumCount, vtkIdType);

  // Description:
  // Specifies a minimum threshold that an edge weight must exceed to be included in
  // the output.
  vtkGetMacro(MinimumThreshold, double);
  vtkSetMacro(MinimumThreshold, double);

protected:
  vtkAdjacencyMatrixToEdgeTable();
  ~vtkAdjacencyMatrixToEdgeTable();

  int FillInputPortInformation(int, vtkInformation*);

  int RequestData(
    vtkInformation*, 
    vtkInformationVector**, 
    vtkInformationVector*);

  vtkIdType SourceDimension;
  char* ValueArrayName;
  vtkIdType MinimumCount;
  double MinimumThreshold;

private:
  vtkAdjacencyMatrixToEdgeTable(const vtkAdjacencyMatrixToEdgeTable&); // Not implemented
  void operator=(const vtkAdjacencyMatrixToEdgeTable&);   // Not implemented
};

#endif

