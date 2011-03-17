/*=========================================================================

 Program:   Visualization Toolkit
 Module:    vtkAMRDualMeshExtractor.h

 Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
 All rights reserved.
 See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

 =========================================================================*/
// .NAME vtkAMRDualMeshExtractor.h -- Extracts the dual mesh from an AMR dataset
//
// .SECTION Description
// This is a concrete instance of a vtkMultiBlockDataSetAlgorithm which accepts
// as input an AMR dataset, represented in a vtkHierarchicalBoxDataSet instance,
// and outputs the dual-mesh of each block given in corresponding instance of
// vtkMultiBlockDataSet.

#ifndef VTKAMRDUALMESHEXTRACTOR_H_
#define VTKAMRDUALMESHEXTRACTOR_H_

#include "vtkMultiBlockDataSetAlgorithm.h"

// Forward declarations
class vtkInformation;
class vtkInformationVector;
class vtkUniformGrid;
class vtkHierarchicalBoxDataSet;
class vtkMultiBlockDataSet;
class vtkUnstructuredGrid;
class vtkIdList;

class VTK_AMR_EXPORT vtkAMRDualMeshExtractor :
  public vtkMultiBlockDataSetAlgorithm
{
  public:
    static vtkAMRDualMeshExtractor *New();
    vtkTypeMacro(vtkAMRDualMeshExtractor,vtkMultiBlockDataSetAlgorithm);
    void PrintSelf( std::ostream& oss, vtkIndent indent);

    // Description:
    // This method writes multiblock data. Note, this method is mostly used
    // for debugging purposes.
    void WriteMultiBlockData( vtkMultiBlockDataSet *mbds );

  protected:
    vtkAMRDualMeshExtractor();
    virtual ~vtkAMRDualMeshExtractor();

    // Description:
    // TODO: enter description
    void GetNeighbor(
        const int ijk[3], const int dims[3],
        const int di, const int dj, const int dk,
        vtkIdList *neiList );

    // Description:
    // TODO: enter descripion
    void GetCellNeighbors(
       const int cellijk[3], const int celldims[3],
       vtkIdList *neisIdList );

    // Description:
    // TODO: enter description
    bool ProcessCellDual(
     vtkUniformGrid *ug, const int cellIdx,
     const int cellijk[3], const int celldims[3] );

    // Description:
    // This methnod computes the cell point ids for a given ijk point.
    // It Returns true if a valid cell can be formed from the given point
    // else false. Note, a valid cell cannot be formed if the point is on
    // a max boundary w.r.t. to the given dimensions.
    bool GetCellIds(  int ijk[3], int dims[3],
          vtkIdList *pntIdList, int numNodesPerCell );

    // Description:
    // This method computes the center of the given cell.
    void ComputeCellCenter(vtkUniformGrid *ug, int cellIdx, double center[3]);

    // Description:
    // TODO: enter description
    void ExtractDualMesh(
     vtkHierarchicalBoxDataSet *amrds, vtkMultiBlockDataSet* mbds );

    // Descrition:
    // TODO: enter description
    vtkUnstructuredGrid* GetDualMesh( vtkUniformGrid *ug );

    // Standard pipeline routines
    virtual int RequestData(
        vtkInformation*,vtkInformationVector**,vtkInformationVector*);
    virtual int FillInputPortInformation(int port, vtkInformation *info);
    virtual int FillOutputPortInformation(int port, vtkInformation *info);

  private:
    vtkAMRDualMeshExtractor(const vtkAMRDualMeshExtractor&); // Not implemented
    void operator=(const vtkAMRDualMeshExtractor&); // Not implemented
};

#endif /* VTKAMRDUALMESHEXTRACTOR_H_ */
