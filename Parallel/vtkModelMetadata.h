/*=========================================================================

  Program:   ParaView
  Module:    vtkModelMetadata.h

  Copyright (c) Kitware, Inc.
  All rights reserved.
  See Copyright.txt or http://www.paraview.org/HTML/Copyright.html for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
/*----------------------------------------------------------------------------
 Copyright (c) Sandia Corporation
 See Copyright.txt or http://www.paraview.org/HTML/Copyright.html for details.
----------------------------------------------------------------------------*/

// .NAME vtkModelMetadata
//
// .SECTION Description
//   This class encapsulates the initialization and
//   static model data that is in an Exodus II file but which is not
//   stored in a vtkUnstructuredGrid.  It can be initialized by
//   the Exodus reader (using the vtkExodusModel class, which
//   interacts with the Exodus library).  It can be supplied 
//   later on to the Exodus writer to improve the quality of the output.  
//
//   Because this class does not depend on the Exodus library, it 
//   should be possible to use it to represent metadata for other 
//   dataset file formats.
//  
//   The fields in this class are those described in the document
//   "EXODUS II: A Finite Element Data Model", SAND92-2137, November 1995.
//
//   Element and node IDs stored in this object must be global IDs.  
//   (The IDs stored in the Exodus file are "internal IDs", or indices into 
//   their location in an array in the Exodus file.  They will be
//   replaced with internal IDs if the file is written back out.)
//
//   One way to initialize this object is by using vtkExodusModel.
//   That class will take an open Exodus II file and a
//   vtkUnstructuredGrid drawn from it and will set the required fields.
//
//   Alternatively, you can use all the Set*
//   methods to set the individual fields. This class does not 
//   copy the data, it simply uses your pointer. This 
//   class will free the storage associated with your pointer 
//   when the class is deleted.  Most fields have sensible defaults.
//   The only requirement is that if you are using this ModelMetadata
//   to write out an Exodus file in parallel, you must SetBlockIds and
//   SetBlockIdArrayName.  Your vtkUnstructuredGrid must have a
//   cell array giving the block ID for each cell.
//
// .SECTION Caveats
//   It is likely that the vtkUnstructuredGrid created from the Exodus
//   file or files will be subsetted, redistributed, or perhaps will
//   aquire ghost cells, prior to being written out by the writer.
//   So it is essential to have global node and element IDs, and a 
//   block ID array in the output grid, so these can be mapped to 
//   tables in this model. 
//
//   The Exodus II library supports an optimized element order map 
//   (section 3.7 in the SAND document).  It contains all the element 
//   IDs, listed in the order in which a solver should process them.  
//   We don't include this, and won't unless there is a request.
//
// .SECTION See also
//   vtkExodusModel vtkExodusReader  vtkExodusIIWriter 

#ifndef __vtkModelMetadata_h
#define __vtkModelMetadata_h

#include <vtkstd/set>
#include <vtkstd/map>

#include "vtkObject.h"

#define myVtkGetMacro(name, type) virtual type Get##name() const { return this->name; }

#define myVtkGetStringMacro(name) virtual char* Get##name () const { return this->name; }

class vtkDataSet;
class vtkCharArray;
class vtkIntArray;
class vtkFloatArray;
class vtkIntArray;

class VTK_EXPORT vtkModelMetadata : public vtkObject
{ 
public:
  vtkTypeRevisionMacro(vtkModelMetadata, vtkObject);
  virtual void PrintSelf(ostream &os, vtkIndent indent);
  static vtkModelMetadata *New();

  // Description:
  //    The global fields are those which pertain to the whole
  //    file.  Examples are the title, information lines,
  //    and list of block IDs.  This method prints out all the
  //    global information.

  virtual void PrintGlobalInformation();

  // Description:
  //    The local fields are those which depend on exactly which
  //    blocks, which time step, and which variables you read in
  //    from the file.  Examples are the number of cells in
  //    each block, and the list of nodes in a node set, or the
  //    value of the global variables at a time step.  If
  //    VERBOSE_TESTING is defined in your execution environment,
  //    this method will print more than mere counts, and actually
  //    print a few of the IDs, distribution factors and so on.  If
  //    VERY_VERBOSE_TESTING is defined, it will print out 
  //    all ID lists, distribution factor lists, and so on.

  virtual void PrintLocalInformation();

  // Description:
  //   The title of the dataset.
  vtkSetStringMacro(Title);
  myVtkGetStringMacro(Title);

  // Description:
  //   Set the information lines. At most 
  //   MAX_LINE_LENGTH characters in each line will be used.
  void SetInformationLines(int numLines, char **lines);
  
  // Description:
  //   Add an information line.  The number 
  //   of characters written out to the file will not exceed
  //   MAX_LINE_LENGTH.  
  void AddInformationLine(char *info);

  // Description:
  //   Get a pointer to all the information lines.  The number
  //   of lines is returned;
  int GetInformationLines(char ***lines) const;

  // Description:
  //   Get the number of information lines.
  myVtkGetMacro(NumberOfInformationLines, int);

  // Description:
  //   Set the list of QA records.  If there was already a
  //   a list, it will be replaced with this one.  We use your
  //   pointer and delete the records when done.
  void SetQARecords(int numberOfRecords, char *QARecords[][4]); 
  
  // Description:
  //   Add a QA record.  Each of the 4 fields must no larger
  //   that MAX_STR_LENGTH.  They are:
  //    The code name
  //    The code version number
  //    The date (MM/DD/YY or NULL for today)
  //    The time (HH:MM:SS or NULL for right now)
  void AddQARecord(char *name, char *version, char *date, char *time);

  // Description:
  //   Get a pointer to the 4 fields of a QA record
  void GetQARecord(int which, 
          char **name, char **version, char **date, char **time) const;

  // Description:
  //   Get the number of QA records
  myVtkGetMacro(NumberOfQARecords, int);

  // Description:
  //    Set the index of the time step represented by the results
  //    data in the file attached to this ModelMetadata object.  Time
  //    step indices start at 0 in this file, they start at 1 in
  //    an Exodus file.
  vtkSetMacro(TimeStepIndex, int);
  myVtkGetMacro(TimeStepIndex, int);

  // Description:
  //    Set the total number of time steps in the file,
  //    and the value at each time step.  We use your time
  //    step value array and delete it when we're done.
  void SetTimeSteps(int numberOfTimeSteps, float *timeStepValues);
  myVtkGetMacro(NumberOfTimeSteps, int);

  // Description:
  //    Get the time step values
  float *GetTimeStepValues() const {return this->TimeStepValues;}

  // Description:
  //   The name of the one, two or three coordinate dimensions.
  //   Each name saved to the file will not exceed MAX_STR_LENGTH
  //   characters.
  void SetCoordinateNames(int dimension, char **);
  char **GetCoordinateNames() const {return this->CoordinateNames;}

  // Description:
  //   Get the dimension of the model.  This is also the number
  //   of coordinate names.
  myVtkGetMacro(Dimension, int);

  // Description:
  //   The number of blocks in the file.  Set this before setting
  //   any of the block arrays.
  vtkSetMacro(NumberOfBlocks, int);
  myVtkGetMacro(NumberOfBlocks, int);

  // Description:
  //   An arbitrary integer ID for each block.
  //   We use your pointer, and free the memory when the object is freed.
  void SetBlockIds(int *);
  int *GetBlockIds() const {return this->BlockIds;}

  // Description:
  //   Element type for each block - a name that means 
  //   something to person who created the file.  At most 
  //   MAX_STR_LENGTH characters per name are written to file.
  //   We use your pointers, and free the memory when the object is freed.
  void SetBlockElementType(char **);
  char **GetBlockElementType() const {return this->BlockElementType;}

  // Description:
  //   Set or get a pointer to a list of the number of elements in
  //   each block.
  //   We use your pointers, and free the memory when the object is freed.
  int SetBlockNumberOfElements(int *nelts);
  int *GetBlockNumberOfElements()const{return this->BlockNumberOfElements;}

  // Description:
  //   Set or get a pointer to a list of the number of nodes in the
  //   elements of  each block.
  //   We use your pointers, and free the memory when the object is freed.
  void SetBlockNodesPerElement(int *);
  int *GetBlockNodesPerElement()const{return this->BlockNodesPerElement;}

  // Description:
  //   Set or get a pointer to a list global element IDs for the
  //   elements in each block. 
  //   We use your pointers, and free the memory when the object is freed.
  void SetBlockElementIdList(int *);
  int *GetBlockElementIdList() const {return this->BlockElementIdList;}

  // Description:
  //    Get the length of the list of elements in every block. 
  myVtkGetMacro(SumElementsPerBlock, int);                                         

  // Description:
  //   Get a list of the index into the BlockElementIdList of the
  //   start of each block's elements.
  int *GetBlockElementIdListIndex()const {return this->BlockElementIdListIndex;} 

  // Description:
  //   Set or get a pointer to a list of the number of attributes
  //   stored for the elements in each block.
  //   We use your pointers, and free the memory when the object is freed.
  int SetBlockNumberOfAttributesPerElement(int *natts);
  int *GetBlockNumberOfAttributesPerElement()const {return this->BlockNumberOfAttributesPerElement;}

  // Description:
  //    Set or get a pointer to a list of the attributes for all
  //    blocks.  The order of the list should be by block, by element
  //    within the block, by attribute.  Omit blocks that don't
  //    have element attributes.
  void SetBlockAttributes(float *);
  float *GetBlockAttributes()const {return this->BlockAttributes;}

  // Description:
  //    Get the length of the list of floating point block attributes.
  myVtkGetMacro(SizeBlockAttributeArray, int);

  // Description:
  //   Get a list of the index into the BlockAttributes of the
  //   start of each block's element attribute list.
  int *GetBlockAttributesIndex()const {return this->BlockAttributesIndex;};

  // Description:
  //   The number of node sets in the file.  Set this value before
  //   setting the various node set arrays.
  vtkSetMacro(NumberOfNodeSets, int);
  myVtkGetMacro(NumberOfNodeSets, int);

  // Description:
  //   Set or get the list the IDs for each node set.
  //   Length of list is the number of node sets.
  //   We use your pointer, and free the memory when the object is freed.
  void SetNodeSetIds(int *);
  int *GetNodeSetIds()const {return this->NodeSetIds;}

  // Description:
  //   Set or get a pointer to a list of the number of nodes in each node set.
  //   We use your pointer, and free the memory when the object is freed.
  int SetNodeSetSize(int *);
  int *GetNodeSetSize()const {return this->NodeSetSize;}

  // Description:
  //   Set or get a pointer to a concatenated list of the
  //   IDs of all nodes in each node set.  First list all IDs in
  //   node set 0, then all IDs in node set 1, and so on.
  //   We use your pointer, and free the memory when the object is freed.
  void SetNodeSetNodeIdList(int *);
  int *GetNodeSetNodeIdList()const {return this->NodeSetNodeIdList;}

  // Description:
  //   Set or get a list of the number of distribution factors stored
  //   by each node set.  This is either 0 or equal to the number of
  //   nodes in the node set. 
  //   Length of list is number of node sets.
  //   We use your pointer, and free the memory when the object is freed.
  int SetNodeSetNumberOfDistributionFactors(int *);
  int *GetNodeSetNumberOfDistributionFactors()const {return this->NodeSetNumberOfDistributionFactors;}

  // Description:
  //   Set or get a list of the distribution factors for the node sets.
  //   The list is organized by node set, and within node set by node.
  //   We use your pointer, and free the memory when the object is freed.
  void SetNodeSetDistributionFactors(float *);
  float *GetNodeSetDistributionFactors()const {return this->NodeSetDistributionFactors;}

  // Description:
  //   Get the total number of nodes in all node sets
  myVtkGetMacro(SumNodesPerNodeSet, int);     

  // Description:
  //   Get the total number of distribution factors stored for all node sets
  myVtkGetMacro(SumDistFactPerNodeSet, int);

  // Description:
  //   Get a list of the index of the starting entry for each node set
  //   in the list of node set node IDs.
  int *GetNodeSetNodeIdListIndex() const {return this->NodeSetNodeIdListIndex;}

  // Description:
  //   Get a list of the index of the starting entry for each node set
  //   in the list of node set distribution factors.
  int *GetNodeSetDistributionFactorIndex() const {return this->NodeSetDistributionFactorIndex;}

  // Description:
  //   Set or get the number of side sets.  Set this value before
  //   setting any of the other side set arrays.
  vtkSetMacro(NumberOfSideSets, int);
  myVtkGetMacro(NumberOfSideSets, int);

  // Description:
  //   Set or get a pointer to a list giving the ID of each side set. 
  //   We use your pointer, and free the memory when the object is freed.
  void SetSideSetIds(int *);
  int *GetSideSetIds()const {return this->SideSetIds;}

  // Description:
  //   Set or get a pointer to a list of the number of sides  in each side set.
  //   We use your pointer, and free the memory when the object is freed.
  int SetSideSetSize(int *sizes);
  int *GetSideSetSize()const {return this->SideSetSize;}

  // Description:
  //   Set or get a pointer to a list of the number of distribution
  //   factors stored by each side set.   Each side set has either
  //   no distribution factors, or 1 per node in the side set.
  //   We use your pointer, and free the memory when the object is freed.
  int SetSideSetNumberOfDistributionFactors(int *df);
  int *GetSideSetNumberOfDistributionFactors()const {return this->SideSetNumberOfDistributionFactors;}

  // Description:
  //   Set or get a pointer to a list of the elements containing each
  //   side in each side set.  The list is organized by side set, and
  //   within side set by element, by node. IS THAT RIGHT????
  //   We use your pointer, and free the memory when the object is freed.
  void SetSideSetElementList(int *);
  int *GetSideSetElementList()const {return this->SideSetElementList;}

  // Description:
  //   Set or get a pointer to the element side for each side in the side set.
  //   (See the manual for the convention for numbering sides in different
  //   types of cells.)  Side Ids are arranged by side set and within
  //   side set by side, and correspond to the SideSetElementList.
  //   We use your pointer, and free the memory when the object is freed.
  void SetSideSetSideList( int *);
  int *GetSideSetSideList()const {return this->SideSetSideList;}

  // Description:
  //   Set or get a pointer to a list of the number of nodes in each
  //   side of each side set.  This list is organized by side set, and
  //   within side set by side.
  //   We use your pointer, and free the memory when the object is freed.
  void SetSideSetNumDFPerSide(int *numNodes);
  int *GetSideSetNumDFPerSide()const {return this->SideSetNumDFPerSide;}

  // Description:
  //   Set or get a pointer to a list of all the distribution factors.
  //   For every side set that has distribution factors, the number of
  //   factors per node was given in the SideSetNumberOfDistributionFactors
  //   array.  If this number for a given side set is N, then for that
  //   side set we have N floating point values for each node for each
  //   side in the side set.  If nodes are repeated in more than one
  //   side, we repeat the distribution factors.  So this list is in order
  //   by side set, by node.
  //   We use your pointer, and free the memory when the object is freed.
  void SetSideSetDistributionFactors(float *);
  float *GetSideSetDistributionFactors()const {return this->SideSetDistributionFactors;}

  // Description:
  //   Get the total number of sides in all side sets
  myVtkGetMacro(SumSidesPerSideSet, int);

  // Description:
  //   Get the total number of distribution factors stored for all side sets
  myVtkGetMacro(SumDistFactPerSideSet, int);

  // Description:
  //   Get a list of the index of the starting entry for each side set
  //   in the list of side set side IDs.
  int *GetSideSetListIndex()const {return this->SideSetListIndex;}

  // Description:
  //   Get a list of the index of the starting entry for each side set
  //   in the list of side set distribution factors.
  int *GetSideSetDistributionFactorIndex()const {return this->SideSetDistributionFactorIndex;}

  // Description:
  //   The number of block properties (global variables)
  myVtkGetMacro(NumberOfBlockProperties, int);

  // Description:
  //   Set or get the names of the block properties.  At most
  //   MAX_STR_LENGTH characters are saved for each name.
  void SetBlockPropertyNames(int numProp, char **names);
  char **GetBlockPropertyNames()const {return this->BlockPropertyNames;}

  // Description:
  //   Set or get value for each variable for each block.  List
  //   the integer values in order by variable and within variable
  //   by block.
  void SetBlockPropertyValue(int *);
  int *GetBlockPropertyValue()const {return this->BlockPropertyValue;}

  // Description:
  //   The number of node set properties (global variables)
  myVtkGetMacro(NumberOfNodeSetProperties, int);

  // Description:
  //   Set or get the names of the node setproperties.  At most
  //   MAX_STR_LENGTH characters are saved for each name.
  void SetNodeSetPropertyNames(int numProp, char **names);
  char **GetNodeSetPropertyNames()const {return this->NodeSetPropertyNames;}

  // Description:
  //   Set or get value for each variable for each node set.  List
  //   the integer values in order by variable and within variable
  //   by node set.
  void SetNodeSetPropertyValue(int *);
  int *GetNodeSetPropertyValue()const {return this->NodeSetPropertyValue;}

  // Description:
  //   The number of side set properties (global variables)
  myVtkGetMacro(NumberOfSideSetProperties, int);

  // Description:
  //   Set or get the names of the side set properties.  At most
  //   MAX_STR_LENGTH characters are saved for each name.
  void SetSideSetPropertyNames(int numProp, char **names);
  char **GetSideSetPropertyNames()const {return this->SideSetPropertyNames;}

  // Description:
  //   Set or get value for each variable for each side set.  List
  //   the integer values in order by variable and within variable
  //   by side set.
  void SetSideSetPropertyValue(int *);
  int *GetSideSetPropertyValue()const {return this->SideSetPropertyValue;}

  // Description:
  //   Get the number of global variables per time step
  myVtkGetMacro(NumberOfGlobalVariables, int);

  // Description:
  //  Set or get the names of the global variables
  void SetGlobalVariableNames(int numVarNames, char **n);
  char **GetGlobalVariableNames()const {return this->GlobalVariableNames;}

  // Description:
  //   Set or get the values of the global variables at the current
  //   time step.
  void SetGlobalVariableValue(float *f);
  float *GetGlobalVariableValue()const {return this->GlobalVariableValue;}

  // Description:
  //   The ModelMetadata maintains a list of the element variables that
  //   were in the original file, and a list of the cell variables
  //   in the UGrid derived from that file.  Some of the scalar variables
  //   in the original file were combined into vectors in the UGrid.
  //   In this method, provide the number of original element variables,
  //   the names of the original element variables, the number of
  //   element variables in the UGrid, the number of components for each
  //   of those variables, and a map from each UGrid variable to the
  //   the variable in the list of original names that represents it's
  //   first component.
  void SetElementVariableInfo(int numOrigNames, char **origNames,
            int numNames, char **names,  int *numComp, int *map);

  // Description:
  //   The ModelMetadata maintains a list of the node variables that
  //   were in the original file, and a list of the node variables
  //   in the UGrid derived from that file.  Some of the scalar variables
  //   in the original file were combined into vectors in the UGrid.
  //   In this method, provide the number of original node variables,
  //   the names of the original node variables, the number of
  //   node variables in the UGrid, the number of components for each
  //   of those variables, and a map from each UGrid variable to the
  //   the variable in the list of original names that represents it's
  //   first component.
  void SetNodeVariableInfo(int numOrigNames, char **origNames,
            int numNames, char **names,  int *numComp, int *map);

  // Description:
  //   A truth table indicating which element variables are
  //   defined for which blocks. The variables are all the original
  //   element variables that were in the file.
  //   The table is by block ID and within block ID by variable.
  void SetElementVariableTruthTable(int *);
  int *GetElementVariableTruthTable()const {return this->ElementVariableTruthTable;}

  // Description:
  //   Instead of a truth table of all "1"s, you can set this
  //   instance variable to indicate that all variables are
  //   defined in all blocks.
  vtkSetMacro(AllVariablesDefinedInAllBlocks, int);
  myVtkGetMacro(AllVariablesDefinedInAllBlocks, int);
  vtkBooleanMacro(AllVariablesDefinedInAllBlocks, int);

  // Description:
  //   If the element variable named is defined for the block Id
  //   provided (in the element variable truth table) return a
  //   1, otherwise return a 0.  If the variable name or block Id
  //   are unrecognized, the default value of 1 is returned.
  //   (This is an "original" variable name, from the file,
  //   not a name created for the vtkUnstructuredGrid.  Use
  //   FindOriginal*VariableName to map between the two.)
  int ElementVariableIsDefinedInBlock(char *varname, int blockId);

  // Description:
  //   The ModelMetadata object may contain these lists:
  //    o  the variables in the original data file
  //    o  the variables created in the u grid from those original variables
  //    o  a mapping from the grid variable names to the original names
  //    o  a list of the number of components each grid variable has
  //
  //   (Example: Variables in Exodus II files are all scalars.  Some are
  //   combined by the ExodusReader into vector variables in the grid.)
  //
  //   These methods return names of the original variables, the names
  //   of the grid variables, a list of the number of components in
  //   each grid variable, and a list of the index into the list of
  //   original variable names where the original name of the first 
  //   component of a grid variable may be found.  The names of subsequent
  //   components would immediately follow the name of the the first
  //   component.
  myVtkGetMacro(OriginalNumberOfElementVariables, int);
  char **GetOriginalElementVariableNames()const {return this->OriginalElementVariableNames;}
  myVtkGetMacro(NumberOfElementVariables, int);
  char **GetElementVariableNames()const {return this->ElementVariableNames;}
  int *GetElementVariableNumberOfComponents()const {return this->ElementVariableNumberOfComponents;}
  int *GetMapToOriginalElementVariableNames()const {return this->MapToOriginalElementVariableNames;}

  myVtkGetMacro(OriginalNumberOfNodeVariables, int);
  char **GetOriginalNodeVariableNames()const {return this->OriginalNodeVariableNames;}
  myVtkGetMacro(NumberOfNodeVariables, int);
  char **GetNodeVariableNames()const {return this->NodeVariableNames;}
  int *GetNodeVariableNumberOfComponents()const {return this->NodeVariableNumberOfComponents;}
  int *GetMapToOriginalNodeVariableNames()const {return this->MapToOriginalNodeVariableNames;}

  // Description:
  //   Given the name of an element variable the vtkUnstructuredGrid
  //   described by this ModelMetadata, and a component number, give 
  //   the name of the scalar array in the original
  //   file that turned into that component when the file was
  //   read into VTK.
  char *FindOriginalElementVariableName(const char *name, int component);

  // Description:
  //   Given the name of an node variable the vtkUnstructuredGrid
  //   described by this ModelMetadata, and a component number, give 
  //   the name of the scalar array in the original
  //   file that turned into that component when the file was
  //   read into VTK.
  char *FindOriginalNodeVariableName(const char *name, int component);

  // Description:
  //   Static function that returns 1 if the vtkUnstructuredGrid
  //   has metadata packed into it's field arrays, and 0 otherwise.
  static int HasMetadata(vtkDataSet *grid);

  // Description:
  //   Pack this object's metadata into a field array of a dataset.
  void Pack(vtkDataSet *ugrid);

  // Description:
  //   Unpack the metadata stored in a dataset,
  //   and initialize this object with it.  Return 1 if there's
  //   no metadata packed into the grid, 0 if OK.
  //   If deleteIt is ON, then delete the grid's packed data after
  //   unpacking it into the object.
  int Unpack(vtkDataSet *ugrid, int deleteIt);

  // Description:
  //   In order to write Exodus files from vtkUnstructuredGrid
  //   objects that were read from Exodus files, we need to know
  //   the mapping from variable names in the UGrid to variable
  //   names in the Exodus file.  (The Exodus reader combines
  //   scalar variables with similar names into vectors in the
  //   UGrid.)  When building the UGrid to which this
  //   ModelMetadata refers, add each element and node variable
  //   name with this call, including the name of original variable
  //   that yielded it's first component, and the number of components.
  //   If a variable is removed from the UGrid, remove it from
  //   the ModelMetadata.  (If this information is missing or
  //   incomplete, the ExodusIIWriter can still do something
  //   sensible in creating names for variables.)
  int AddUGridElementVariable(char *ugridVarName, char *origName, int numComponents);
  int RemoveUGridElementVariable(char *ugridVarName);

  int AddUGridNodeVariable(char *ugridVarName, char *origName, int numComponents);
  int RemoveUGridNodeVariable(char *ugridVarName);

  // Description:
  //   In VTK we take vtkUnstructuredGrids and perform
  //   operations on them, including subsetting and merging
  //   grids.  We need to modify the metadata object 
  //   when this happens.  MergeModelMetadata merges the supplied
  //   model (both global and local metadata) into this model.  
  //   The models must be from the same file set.
  //
  //   MergeModelMetadata assumes that no element in one metadata
  //   object appears in the other.  (It doesn't test for duplicate
  //   elements when merging the two metadata objects.) 
  int MergeModelMetadata(const vtkModelMetadata *em);

  // Description:
  //   The metadata is divided into global metadata and local
  //   metadata.  MergeGlobalInformation merges just the
  //   global metadata of the supplied object into the
  //   global metadata of this object.
  int MergeGlobalInformation(const vtkModelMetadata *em);

  // Description:
  //   Create and return a new metadata object which contains
  //   the information for the subset of global cell IDs provided.
  //   We need the grid containing the cells so we can find point
  //   Ids as well, and also the name of the global cell ID array
  //   and the name of the global point ID array.
  vtkModelMetadata *ExtractModelMetadata(vtkIntArray *globalCellIdList,
                                     vtkDataSet *grid,
                                     const char *globalCellIdArrayName,
                                     const char *globalNodeIdArrayName);

  // Description:
  //   Create and return a new metadata object containing only the
  //   global metadata of this metadata object.  
  vtkModelMetadata *ExtractGlobalMetadata();

  // Description:
  //   Free selected portions of the metadata when updating values
  //   in the vtkModelMetadata object.  Resetting a particular field,
  //   (i.e. SetNodeSetIds) frees the previous setting, but if you
  //   are not setting every field, you may want to do a wholesale
  //   "Free" first.
  //
  //   FreeAllGlobalData frees all the fields which don't depend on
  //     which time step, which blocks, or which variables are in the input.
  //   FreeAllLocalData frees all the fields which do depend on which
  //     time step, blocks or variables are in the input.
  //   FreeBlockDependentData frees all metadata fields which depend on
  //     which blocks were read in.
  void FreeAllGlobalData();
  void FreeAllLocalData();
  void FreeBlockDependentData();
  void FreeOriginalElementVariableNames();
  void FreeOriginalNodeVariableNames();
  void FreeUsedElementVariableNames();
  void FreeUsedNodeVariableNames();
  void FreeUsedElementVariables();
  void FreeUsedNodeVariables();

  // Description:
  //   Set the object back to it's initial state
  void Reset();

  // Description:
  //   Block information is stored in arrays.  This method returns
  //   the array index for a given block ID.
  int GetBlockLocalIndex(int id);

protected:

  vtkModelMetadata();
  ~vtkModelMetadata();

private:

  void InitializeAllMetadata();
  void InitializeAllIvars();

  void FreeAllMetadata();
  void FreeAllIvars();

  void FreeQARecords();

  int BuildBlockElementIdListIndex();
  int BuildBlockAttributesIndex();
  int BuildNodeSetNodeIdListIndex();
  int BuildNodeSetDistributionFactorIndex();
  int BuildSideSetListIndex();
  int BuildSideSetDistributionFactorIndex();

  int InitializeFromSizeArray(vtkIntArray *ia, int &maxStr, int &maxLine);
  vtkIntArray *PackSizeArray(int maxStr, int maxLine);
  int InitializeFromIntArray(vtkModelMetadata *sizes, vtkIntArray *ia);
  vtkIntArray *PackIntArray();
  int InitializeFromCharArray(vtkModelMetadata *sizes, 
                   vtkCharArray *uca, int maxStr, int maxLine);
  vtkCharArray *PackCharArray(int maxStr, int maxLine);
  int InitializeFromFloatArray(vtkFloatArray *fa);
  vtkFloatArray *PackFloatArray();

  static char *StrDupWithNew(const char *s);

  static char *WriteLines(char *p, int maxLines, int maxLen, char **lines);
  static char *ReadLines(char ***to, int maxLines, 
                            int maxLen, char *from);
  static char **CopyLines(char **lines, int num);
  static int *CopyInts(int *vals, int num);

  static int FindNameOnList(char *name, char **list, int listLen);

  int MergeIdLists(int numSubLists,
    int *id1, int *id1Idx, int id1Len,
      float *dist1, int *dist1Idx, int dist1Len,
    int *id2, int *id2Idx, int id2Len,
      float *dist2, int *dist2Idx, int dist2Len,
    int **idNew, int **idNewIdx, int *idNewLen,
      float **distNew, int **distNewIdx, int *distNewLen);

  int AppendFloatLists(int numSubLists,
    float *id1, int *id1Idx, int id1Len,
    float *id2, int *id2Idx, int id2Len,
    float **idNew, int **idNewIdx, int *idNewLen);

  int AppendIntegerLists(int numSubLists,
    int *id1, int *id1Idx, int id1Len,
    int *id2, int *id2Idx, int id2Len,
    int **idNew, int **idNewIdx, int *idNewLen);

//BTX
  void ExtractCellsFromBlockData(vtkstd::set<int> *idset, vtkModelMetadata *mmd);
  void ExtractNodesFromNodeSetData(vtkstd::set<int> *idset, vtkModelMetadata *mmd);
  void ExtractSidesFromSideSetData(vtkstd::set<int> *idset, vtkModelMetadata *mmd);
//ETX

  void ShowFloats(char *what, int num, float *f);
  void ShowLines(char *what, int num, char **l);
  void ShowIntArray(char *what, int numx, int numy, int *id);
  void ShowInts(char *what, int num, int *id);
  void ShowListsOfInts(char *what, int *list,
                       int nlists, int *idx, int len, int verbose);
  void ShowListsOfFloats(char *what, float *list,
                         int nlists, int *idx, int len, int verbose);

  void SetOriginalElementVariableNames(int nvars, char **names);
  void SetElementVariableNames(int nvars, char **names);
  void SetElementVariableNumberOfComponents(int *comp);
  void SetMapToOriginalElementVariableNames(int *map);

  void SetOriginalNodeVariableNames(int nvars, char **names);
  void SetNodeVariableNames(int nvars, char **names);
  void SetNodeVariableNumberOfComponents(int *comp);
  void SetMapToOriginalNodeVariableNames(int *map);

  int CalculateMaximumLengths(int &maxString, int &maxLine);

  // Fields in Exodus II file and their size (defined in exodusII.h)
  //   (G - global fields, the same in every file)
  //   (L - local fields, they differ depending on which cells and nodes are
  //        in a file)

  char *Title;                 // MAX_LINE_LENGTH    (G)

  int NumberOfQARecords;       // (G)
//BTX
  char *(*QARecord)[4];        // NumberOfQARecords * 4 * MAX_STR_LENGTH (G)
//ETX

  int NumberOfInformationLines; // (G)
  char **InformationLine;       // each record is MAX_LINE_LENGTH (G)

  int Dimension;            // (G)
  char **CoordinateNames;   // MAX_STR_LENGTH each (at most 3 of these) (G)

  // Time steps

  int TimeStepIndex;     // starting at 0 (Exodus file starts at 1) 
  int NumberOfTimeSteps; // (G)
  float *TimeStepValues; // (G)

  // Block information - arrays that are input with Set*

  int NumberOfBlocks;       // (G)

  int *BlockIds;               // NumberOfBlocks (G) (start at 1)
  char **BlockElementType;     // NumberOfBlocks, length MAX_STR_LENGTH (G)
  int *BlockNumberOfElements;  // NumberOfBlocks (L)
  int *BlockNodesPerElement;   // NumberOfBlocks (G)
  int *BlockNumberOfAttributesPerElement;// NumberOfBlocks (G)
  int *BlockElementIdList;     // SumElementsPerBlock     (L)
  float *BlockAttributes;      // SizeBlockAttributeArray (L)

  // Block information - values that we calculate

  int SumElementsPerBlock;                                         
  int SizeBlockAttributeArray;

  int *BlockElementIdListIndex;          // NumberOfBlocks         
  int *BlockAttributesIndex;             // NumberOfBlocks

//BTX
  vtkstd::map<int, int> BlockIdIndex;
//ETX

  // Node Sets - arrays that are input to the class with Set*

  int NumberOfNodeSets; // (G)

  int *NodeSetIds;             // NumberOfNodeSets (G)
  int *NodeSetSize;            // NumberOfNodeSets (L)
  int *NodeSetNumberOfDistributionFactors;  // NNS (L) (NSNDF[i] is 0 or NSS[i])
  int *NodeSetNodeIdList;   // SumNodesPerNodeSet (L)
  float *NodeSetDistributionFactors; // SumDistFactPerNodeSet (L)

  // Node Sets - values or arrays that the class computes

  int SumNodesPerNodeSet;
  int SumDistFactPerNodeSet;

  int *NodeSetNodeIdListIndex;           // NumberOfNodeSets
  int *NodeSetDistributionFactorIndex;   // NumberOfNodeSets

  // Side Sets - input to class with Set*

  int NumberOfSideSets; // (G)
  
  int *SideSetIds;                          // NumberOfSideSets (G)
  int *SideSetSize;                         // NumberOfSideSets (L)
  int *SideSetNumberOfDistributionFactors;  // NSS (L) (SSNDF[i] = 0 or NumNodesInSide)
  int *SideSetElementList;               // SumSidesPerSideSet (L)
  int *SideSetSideList;                  // SumSidesPerSideSet (L)
  int *SideSetNumDFPerSide;              // SumSidesPerSideSet (L)
  float *SideSetDistributionFactors;     // SumDistFactPerSideSet (L)

  // Side Sets - calculated by class

  int SumSidesPerSideSet;
  int SumDistFactPerSideSet;

  int *SideSetListIndex;                 // NumberOfSideSets
  int *SideSetDistributionFactorIndex;   // NumberOfSideSets

  // Other properties, provided as input with Set*

  int NumberOfBlockProperties; // (G)
  char **BlockPropertyNames;   // one per property, MAX_STR_LENGTH (G)
  int *BlockPropertyValue;     // NumBlocks * NumBlockProperties (G)

  int NumberOfNodeSetProperties; // (G)
  char **NodeSetPropertyNames;   // one per property, MAX_STR_LENGTH (G)
  int *NodeSetPropertyValue;     // NumNodeSets * NumNodeSetProperties (G)

  int NumberOfSideSetProperties; // (G)
  char **SideSetPropertyNames;   // one per property, MAX_STR_LENGTH (G)
  int *SideSetPropertyValue;     // NumSideSets * NumSideSetProperties (G)

  // Global variables, 1 value per time step per variable.  We store
  // these as floats, even if they are doubles in the file.  The values
  // are global in the sense that they apply to the whole data set, but
  // the are local in the sense that they can change with each time step.
  // For the purpose of this object, which represents a particular
  // time step, they are therefore considered "local".  (Since they need
  // to be updated everytime another read is done from the file.)

  int NumberOfGlobalVariables;   // (G)
  char **GlobalVariableNames;    // (G) NumberOfGlobalVariables, MAX_STR_LENGTH
  float *GlobalVariableValue;   // (G) NumberOfGlobalVariables

  // The element and node arrays in the file were all scalar arrays.
  // Those with similar names were combined into vectors in VTK.  Here
  // are all the original names from the Exodus file, the names given
  // the variables in the VTK ugrid, and a mapping from the VTK names
  // to the Exodus names.

  int OriginalNumberOfElementVariables;    // (G)
  char **OriginalElementVariableNames;     // (G) OriginalNumberOfElementVariables
  int NumberOfElementVariables;            // (G)
  int MaxNumberOfElementVariables;         // (G)
  char **ElementVariableNames;             // (G) MaxNumberOfElementVariables
  int *ElementVariableNumberOfComponents;  // (G) MaxNumberOfElementVariables
  int *MapToOriginalElementVariableNames;  // (G) MaxNumberOfElementVariables

  int OriginalNumberOfNodeVariables;       // (G)
  char **OriginalNodeVariableNames;        // (G) OriginalNumberOfNodeVariables
  int NumberOfNodeVariables;               // (G)
  int MaxNumberOfNodeVariables;            // (G)
  char **NodeVariableNames;                // (G) NumberOfNodeVariables
  int *NodeVariableNumberOfComponents;     // (G) NumberOfNodeVariables
  int *MapToOriginalNodeVariableNames;     // (G) NumberOfNodeVariables

  int *ElementVariableTruthTable;  // (G) NumBlocks*OrigNumberOfElementVariables
  int AllVariablesDefinedInAllBlocks;

private:
  vtkModelMetadata(const vtkModelMetadata&); // Not implemented
  void operator=(const vtkModelMetadata&); // Not implemented
};
#endif
