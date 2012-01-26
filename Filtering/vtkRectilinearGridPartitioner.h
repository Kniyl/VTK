/*=========================================================================

 Program:   Visualization Toolkit
 Module:    vtkRectilinearGridPartitioner.h

 Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
 All rights reserved.
 See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

 =========================================================================*/
// .NAME vtkRectilinearGridPartitioner.h -- Partitions a rectilinear grid by RCB
//
// .SECTION Description
//  A concrete implementation of vtkMultiBlockDataSetAlgorithm that provides
//  functionality for partitioning a VTK rectilinear dataset. The partitioning
//  methd used is Recursive Coordinate Bisection (RCB) where each time the
//  longest dimension is split.
//
// .SECTION See Also
//  vtkUniformGridPartitioner vtkStructuredGridPartitioner
#ifndef VTKRECTILINEARGRIDPARTITIONER_H_
#define VTKRECTILINEARGRIDPARTITIONER_H_

#include "vtkMultiBlockDataSetAlgorithm.h"

class vtkInformation;
class vtkInformationVector;
class vtkIndent;

class VTK_FILTERING_EXPORT vtkRectilinearGridPartitioner :
  public vtkMultiBlockDataSetAlgorithm
{
  public:
    static vtkRectilinearGridPartitioner *New();
    vtkTypeMacro(vtkRectilinearGridPartitioner, vtkMultiBlockDataSetAlgorithm );
    void PrintSelf( std::ostream &oss, vtkIndent indent );

    // Description:
    // Set/Get macro for the number of subdivisions.
    vtkGetMacro(NumberOfPartitions,int);
    vtkSetMacro(NumberOfPartitions,int);

    // Description:
    // Set/Get macro for the number of ghost layers.
    vtkGetMacro(NumberOfGhostLayers,int);
    vtkSetMacro(NumberOfGhostLayers,int);

  protected:
    vtkRectilinearGridPartitioner();
    virtual ~vtkRectilinearGridPartitioner();

    // Standard Pipeline methods
    virtual int RequestData(
       vtkInformation*,vtkInformationVector**,vtkInformationVector*);
    virtual int FillInputPortInformation(int port, vtkInformation *info);
    virtual int FillOutputPortInformation(int port, vtkInformation *info);

    int NumberOfPartitions;
    int NumberOfGhostLayers;

  private:
    vtkRectilinearGridPartitioner(const vtkRectilinearGridPartitioner &);
    void operator=(const vtkRectilinearGridPartitioner &);
};

#endif /* VTKRECTILINEARGRIDPARTITIONER_H_ */
