/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkGenericDataSetAlgorithm.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkGenericDataSetAlgorithm - Objects that generate adapted data sets
// .SECTION Description

// vtkGenericDataSetAlgorithm is a convenience class to make writing algorithms
// easier. It is also designed to help transition old algorithms to the new
// pipeline architecture. Ther are some assumptions and defaults made by this
// class you should be aware of. This class defaults such that your filter
// will have one input port and one output port. If that is not the case
// simply change it with SetNumberOfInputPorts etc. See this classes
// constructor for the default. This class also provides a FillInputPortInfo
// method that by default says that all inputs will be GenericDataSet. If that
// isn't the case then please override this method in your subclass. This
// class breaks out the downstream requests into seperate functions such as
// ExecuteData and ExecuteInformation.  For new algorithms you should
// implement RequestData( request, inputVec, outputVec) but for older filters
// there is a default implementation that calls the old ExecuteData(output)
// signature, for even older filters that don;t implement ExecuteData the
// default implementation calls the even older Execute() signature.

#ifndef __vtkGenericDataSetAlgorithm_h
#define __vtkGenericDataSetAlgorithm_h

#include "vtkAlgorithm.h"
#include "vtkGenericDataSet.h" // makes things a bit easier

class vtkDataSet;
class vtkGenericDataSet;

class VTK_FILTERING_EXPORT vtkGenericDataSetAlgorithm : public vtkAlgorithm
{
public:
  vtkTypeRevisionMacro(vtkGenericDataSetAlgorithm,vtkAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Get the output data object for a port on this algorithm.
  vtkGenericDataSet* GetOutput();
  vtkGenericDataSet* GetOutput(int);
  virtual void SetOutput(vtkDataObject* d);

  // Description:
  // see vtkAlgorithm for details
  virtual int ProcessRequest(vtkInformation*,
                             vtkInformationVector**,
                             vtkInformationVector*);

  // this method is not recommended for use, but lots of old style filters
  // use it
  vtkDataObject* GetInput();
  vtkDataObject *GetInput(int port);
  vtkGenericDataSet *GetGenericDataSetInput(int port);

  // Description:
  // Set an input of this algorithm.
  void SetInput(vtkDataObject *);
  void SetInput(int, vtkDataObject*);

  // Description:
  // Add an input of this algorithm.
  void AddInput(vtkDataObject *);
  void AddInput(int, vtkDataObject*);

protected:
  vtkGenericDataSetAlgorithm();
  ~vtkGenericDataSetAlgorithm();

  // convinience method
  virtual int RequestInformation(vtkInformation* request,
                                 vtkInformationVector** inputVector,
                                 vtkInformationVector* outputVector);

  // Description:
  // This is called by the superclass.
  // This is the method you should override.
  // See ProcessRequest for details about arguments and return value.
  virtual int RequestData(vtkInformation* request,
                          vtkInformationVector** inputVector,
                          vtkInformationVector* outputVector);
  
  // Description:
  // This is called by the superclass.
  // This is the method you should override.
  // See ProcessRequest for details about arguments and return value.
  virtual int RequestDataObject(vtkInformation* request,
                                vtkInformationVector** inputVector,
                                vtkInformationVector* outputVector)=0;
  
  // Description:
  // This is called by the superclass.
  // This is the method you should override.
 // See ProcessRequest for details about arguments and return value.
  virtual int RequestUpdateExtent(vtkInformation*,
                                  vtkInformationVector**,
                                  vtkInformationVector*);

  // Description:
  // This method is the old style execute method
  virtual void ExecuteData(vtkDataObject *output);
  virtual void Execute();

  // see algorithm for more info
  virtual int FillOutputPortInformation(int port, vtkInformation* info);
  virtual int FillInputPortInformation(int port, vtkInformation* info);

private:
  vtkGenericDataSetAlgorithm(const vtkGenericDataSetAlgorithm&);  // Not implemented.
  void operator=(const vtkGenericDataSetAlgorithm&);  // Not implemented.
};

#endif
