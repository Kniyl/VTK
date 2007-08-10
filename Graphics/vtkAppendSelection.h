/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkAppendSelection.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkAppendSelection - appends one or more selections together
//
// .SECTION Description
// vtkAppendSelection is a filter that appends one of more selections into
// a single selection.  All selections must have the same content type.

#ifndef __vtkAppendSelection_h
#define __vtkAppendSelection_h

#include "vtkSelectionAlgorithm.h"

class vtkSelection;

class VTK_GRAPHICS_EXPORT vtkAppendSelection : public vtkSelectionAlgorithm
{
public:
  static vtkAppendSelection *New();

  vtkTypeRevisionMacro(vtkAppendSelection,vtkSelectionAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // UserManagedInputs allows the user to set inputs by number instead of
  // using the AddInput/RemoveInput functions. Calls to
  // SetNumberOfInputs/SetInputByNumber should not be mixed with calls
  // to AddInput/RemoveInput. By default, UserManagedInputs is false.
  vtkSetMacro(UserManagedInputs,int);
  vtkGetMacro(UserManagedInputs,int);
  vtkBooleanMacro(UserManagedInputs,int);

  // Description:
  // Add a dataset to the list of data to append. Should not be
  // used when UserManagedInputs is true, use SetInputByNumber instead.
  void AddInput(vtkSelection *);

  // Description:
  // Remove a dataset from the list of data to append. Should not be
  // used when UserManagedInputs is true, use SetInputByNumber (NULL) instead.
  void RemoveInput(vtkSelection *);

  // Description:
  // Get any input of this filter.
//BTX
  vtkSelection *GetInput(int idx);
  vtkSelection *GetInput() { return this->GetInput( 0 ); };
//ETX

  // Description:
  // Directly set(allocate) number of inputs, should only be used
  // when UserManagedInputs is true.
  void SetNumberOfInputs(int num);

  // Set Nth input, should only be used when UserManagedInputs is true.
  void SetInputByNumber(int num, vtkSelection *input);

protected:
  vtkAppendSelection();
  ~vtkAppendSelection();

  // Usual data generation method
  virtual int RequestData(vtkInformation *, 
                          vtkInformationVector **, vtkInformationVector *);
  virtual int FillInputPortInformation(int, vtkInformation *);

 private:
  // hide the superclass' AddInput() from the user and the compiler
  void AddInput(vtkDataObject *)
    { vtkErrorMacro( << "AddInput() must be called with a vtkSelection not a vtkDataObject."); };

  int UserManagedInputs;

private:
  vtkAppendSelection(const vtkAppendSelection&);  // Not implemented.
  void operator=(const vtkAppendSelection&);  // Not implemented.
};

#endif


