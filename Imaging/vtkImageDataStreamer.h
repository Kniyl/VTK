/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkImageDataStreamer.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkImageDataStreamer - Initiates streaming on image data.
// .SECTION Description
// To satisfy a request, this filter calls update on its input
// many times with smaller update extents.  All processing up stream
// streams smaller pieces.

#ifndef __vtkImageDataStreamer_h
#define __vtkImageDataStreamer_h

#include "vtkImageToImageFilter.h"

class vtkExtentTranslator;

class VTK_IMAGING_EXPORT vtkImageDataStreamer : public vtkImageToImageFilter
{
public:
  static vtkImageDataStreamer *New();
  vtkTypeRevisionMacro(vtkImageDataStreamer,vtkImageToImageFilter);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Set how many pieces to divide the input into.
  // void SetNumberOfStreamDivisions(int num);
  // int GetNumberOfStreamDivisions();
  vtkSetMacro(NumberOfStreamDivisions,int);
  vtkGetMacro(NumberOfStreamDivisions,int);
  
  virtual void Update() { this->vtkAlgorithm::Update();};  
  virtual void UpdateWholeExtent() { 
    this->vtkAlgorithm::UpdateWholeExtent();};  
  
  // Description:
  // Get the extent translator that will be used to split the requests
  virtual void SetExtentTranslator(vtkExtentTranslator*);
  vtkGetObjectMacro(ExtentTranslator,vtkExtentTranslator);
  
  // See the vtkAlgorithm for a desciption of what these do
  int ProcessRequest(vtkInformation*,
                     vtkInformationVector**,
                     vtkInformationVector*);

protected:
  vtkImageDataStreamer();
  ~vtkImageDataStreamer();

  // see the vtkAlgorithm for what these methods do
  virtual int FillInputPortInformation(int port, vtkInformation* info);
  virtual int FillOutputPortInformation(int port, vtkInformation* info);

  vtkExtentTranslator *ExtentTranslator;
  int            NumberOfStreamDivisions;
  int            CurrentDivision;
private:
  vtkImageDataStreamer(const vtkImageDataStreamer&);  // Not implemented.
  void operator=(const vtkImageDataStreamer&);  // Not implemented.
};

#endif



