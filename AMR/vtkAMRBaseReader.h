/*=========================================================================

 Program:   Visualization Toolkit
 Module:    vtkAMRBaseReader.h

 Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
 All rights reserved.
 See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

 =========================================================================*/
// .NAME vtkAMRBaseReader.h -- Base class for all AMR Readers
//
// .SECTION Description
// An abstract class that encapsulates common functionality for all AMR readers.

#ifndef VTKAMRBASEREADER_H_
#define VTKAMRBASEREADER_H_

#include "vtkHierarchicalBoxDataSetAlgorithm.h"
#include <vtkstd/vector> // STL vector Header

// Forward Declarations
class vtkHierarchicalBoxDataSet;
class vtkMultiProcessController;
class vtkDataArraySelection;
class vtkCallbackCommand;
class vtkIndent;

class VTK_AMR_EXPORT vtkAMRBaseReader :
  public vtkHierarchicalBoxDataSetAlgorithm
{
  public:
     vtkTypeMacro( vtkAMRBaseReader, vtkHierarchicalBoxDataSetAlgorithm );
     void PrintSelf(std::ostream &os, vtkIndent indent);

     // Description:
     // Initializes the AMR reader.
     // All concrete instances must call this method in their constructor.
     void Initialize();

    // Description:
    // Set/Get a multiprocess-controller for reading in parallel.
    // By default this parameter is set to NULL by the constructor.
    vtkSetMacro( Controller, vtkMultiProcessController* );
    vtkGetMacro( Controller, vtkMultiProcessController* );

    // Description:
    // Set the level, up to which the blocks are loaded.
    vtkSetMacro( MaxLevel,int);

    // Description:
    // Set/Get whether the particles, if any, are loaded
    vtkSetMacro( LoadParticles, int );
    vtkGetMacro( LoadParticles, int );
    vtkBooleanMacro( LoadParticles, int );

    // Description:
    // Get the data array selection tables used to configure which data
    // arrays are loaded by the reader.
    vtkGetObjectMacro(CellDataArraySelection, vtkDataArraySelection);
    vtkGetObjectMacro(PointDataArraySelection, vtkDataArraySelection);

    // Description:
    // Get the number of point or cell arrays available in the input.
    int GetNumberOfPointArrays();
    int GetNumberOfCellArrays();

    // Description:
    // Get the name of the point or cell array with the given index in
    // the input.
    const char* GetPointArrayName(int index);
    const char* GetCellArrayName(int index);

   // Description:
   // Get/Set whether the point or cell array with the given name is to
   // be read.
   int GetPointArrayStatus(const char* name);
   int GetCellArrayStatus(const char* name);
   void SetPointArrayStatus(const char* name, int status);
   void SetCellArrayStatus(const char* name, int status);

    // Description:
    // Set/Get the filename. Concrete instances of this class must implement
    // the SetFileName method accordingly.
    vtkGetStringMacro( FileName );
    virtual void SetFileName( const char *fileName ) = 0;

  protected:
    vtkAMRBaseReader();
    ~vtkAMRBaseReader();

    // Description:
    // Determines if the block is owned by this process based on the
    // the block index and total number of processes.
    bool IsBlockMine( const int blockIdx );

    // Description:
    // Returns the block process ID for the block corresponding to the
    // given block index. If this reader instance is serial, i.e., there
    // is no controller associated, the method returns 0. Otherwise, static
    // block-cyclic-distribution is assumed and each block is assigned to
    // a process according to blockIdx%N, where N is the total number of
    // processes.
    int GetBlockProcessId( const int blockIdx );

    // Description:
    // Reads all the metadata from the file. Implemented by concrete classes.
    virtual void ReadMetaData() = 0;

    // Description:
    // Generates a linear index for each block. Implemented by concrete
    // instances.
    //
    // Note this method must create a block map for all blocks
    // that are to be rendered according to the max level argument. For
    // example, a sample implementation may look as follows:
    //<code>
    //  for(int i=0;i < this->Internal->NumberOfBlocks;i++ )
    //   {
    //    if ( this->GetBlockLevel( i ) <= this->MaxLevel  )
    //     {
    //      this->BlockMap.push_back( i );
    //     }
    //   }
    //</code>
    virtual void GenerateBlockMap() = 0;

    // Description:
    // Returns the block level for the given block
    virtual int GetBlockLevel( const int blockIdx ) = 0;

    // Description:
    // Returns the total number of blocks. Implemented by concrete instances.
    virtual int GetNumberOfBlocks() = 0;

    // Description:
    // Returns the total number of levels. Implemented by concrete instances.
    virtual int GetNumberOfLevels() = 0;

    // Description:
    // Loads the block according to the index w.r.t. the generated BlockMap.
    virtual void GetBlock(
        int index, vtkHierarchicalBoxDataSet *hbds,
        vtkstd::vector< int > &idxcounter ) = 0;

    // Description:
    // Standard Pipeline methods, subclasses may override this method if needed.
    virtual int RequestData(
        vtkInformation* vtkNotUsed(request),
        vtkInformationVector** vtkNotUsed(inputVector),
        vtkInformationVector* outputVector );
    int FillOutputPortInformation(int port,vtkInformation *info);

    // Array selection member variables and methods
    vtkDataArraySelection *PointDataArraySelection;
    vtkDataArraySelection *CellDataArraySelection;
    vtkCallbackCommand    *SelectionObserver;

    // Description:
    // Initializes the PointDataArraySelection & CellDataArraySelection
    virtual void SetUpDataArraySelections() = 0;

    // Descriptions
    // Call-back registered with the SelectionObserver.
    static void SelectionModifiedCallback(
      vtkObject *caller,unsigned long eid,void *clientdata,void *calldata );

    int LoadParticles;
    int MaxLevel;
    char *FileName;
    vtkMultiProcessController *Controller;

    //BTX
      vtkstd::vector<int> BlockMap;
    //ETX

  private:
    vtkAMRBaseReader( const vtkAMRBaseReader& ); // Not implemented
    void operator=( const vtkAMRBaseReader& ); // Not implemented
};

#endif /* VTKAMRBASEREADER_H_ */
