/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkGeoTerrainSource.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
/*-------------------------------------------------------------------------
  Copyright 2008 Sandia Corporation.
  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
  the U.S. Government retains certain rights in this software.
-------------------------------------------------------------------------*/
// .NAME vtkGeoTerrainSource - Supper class for terrain source. 
// .SECTION This should contain the API for terrain sources.

// .SECTION See Also
// vtkGeoTerrainGlobeSource
// vtkGeoTerrainRemoteSource

#ifndef __vtkGeoTerrainSource_h
#define __vtkGeoTerrainSource_h

#include "vtkObject.h"

class vtkGeoTerrainNode;

class VTK_GEOVIS_EXPORT vtkGeoTerrainSource : public vtkObject
{
public:
  static vtkGeoTerrainSource *New();
  vtkTypeRevisionMacro(vtkGeoTerrainSource, vtkObject);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // The poly data is shifted wo this world location is at (0,0,0).
  // This is a fix for the OpenGL imprecision bug seen during picking.
  vtkSetVector3Macro(Origin, double);
  
  // Description:
  // This is used by the local globe source.  It returns
  // when the request has been satisfied.
  // I think for the remote class, we should have a callback that
  // sets the completed node in the tree.
  virtual void GenerateTerrainForNode(vtkGeoTerrainNode*)
    {vtkErrorMacro("Method not implemented"); }
  
protected:
  vtkGeoTerrainSource();
  ~vtkGeoTerrainSource();

  double Origin[3];

private:
  vtkGeoTerrainSource(const vtkGeoTerrainSource&);  // Not implemented.
  void operator=(const vtkGeoTerrainSource&);  // Not implemented.
};

#endif
