/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkDEMReader.h
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

  Copyright (c) 1993-2002 Ken Martin, Will Schroeder, Bill Lorensen 
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkDEMReader - read a digital elevation model (DEM) file
// .SECTION Description
// vtkDEMReader reads digital elevation files and creates image data.
// Digital elevation files are produced by the
// <A HREF="http://www.usgs.org">US Geological Survey</A>. 
// A complete description of the DEM file is located at the USGS site.
// The reader reads the entire dem file and create a vtkImageData that
// contains a single scalar component that is the elevation in meters.
// The spacing is also expressed in meters. A number of get methods
// provide access to fields on the header.
#ifndef __vtkDEMReader_h
#define __vtkDEMReader_h

#include <stdio.h>
#include "vtkImageSource.h"
class VTK_IO_EXPORT vtkDEMReader : public vtkImageSource
{
public:
  static vtkDEMReader *New();
  vtkTypeRevisionMacro(vtkDEMReader,vtkImageSource);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Specify file name of Digital Elevation Model (DEM) file
  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);

  // Description:
  // An ASCII description of the map
  vtkGetStringMacro(MapLabel);

  // Description:
  // Code 1=DEM-1, 2=DEM_2, ...
  vtkGetMacro(DEMLevel,int);

  // Description:
  // Code 1=regular, 2=random, reserved for future use
  vtkGetMacro(ElevationPattern,  int);

  // Description:
  // Ground planimetric reference system
  vtkGetMacro(GroundSystem,  int);

  // Description:
  // Zone in ground planimetric reference system
  vtkGetMacro(GroundZone,  int);

  // Description:
  // Map Projection parameters. All are zero.
  vtkGetVectorMacro(ProjectionParameters,float,15);

  // Description:
  // Defining unit of measure for ground planimetric coordinates throughout
  // the file. 0 = radians, 1 = feet, 2 = meters, 3 = arc-seconds.
  vtkGetMacro(PlaneUnitOfMeasure,  int);

  // Description:
  // Defining unit of measure for elevation coordinates throughout
  // the file. 1 = feet, 2 = meters
  vtkGetMacro(ElevationUnitOfMeasure,  int);

  // Description:
  // Number of sides in the polygon which defines the coverage of
  // the DEM file. Set to 4.
  vtkGetMacro(PolygonSize,  int);

  // Description:
  // Minimum and maximum elevation for the DEM. The units in the file
  // are in ElevationUnitOfMeasure. This class converts them to meters.
  vtkGetVectorMacro(ElevationBounds,float,2);

  // Description:
  // Counterclockwise angle (in radians) from the primary axis of the planimetric
  // reference to the primary axis of the DEM local reference system.
  // IGNORED BY THIS IMPLEMENTATION.
  vtkGetMacro(LocalRotation,  float);

  // Description:
  // Accuracy code for elevations. 0=unknown accuracy
  vtkGetMacro(AccuracyCode,  int);

  // Description:
  // DEM spatial resolution for x,y,z. Values are expressed in units of resolution.
  // Since elevations are read as integers, this permits fractional elevations.
  vtkGetVectorMacro(SpatialResolution,float,3);

  // Description:
  // The number of rows and columns in the DEM.
  vtkGetVectorMacro(ProfileDimension,int,2);

  // Description:
  // Reads the DEM Type A record to compute the extent, origin and
  // spacing of the image data. The number of scalar components is set
  // to 1 and the output scalar type is VTK_FLOAT. 
  void ExecuteInformation();

protected:
  vtkDEMReader();
  ~vtkDEMReader();

  vtkTimeStamp ReadHeaderTime;
  int NumberOfColumns;
  int NumberOfRows;
  int WholeExtent[6];
  char *FileName;
  char MapLabel[145];
  int DEMLevel;
  int ElevationPattern;
  int GroundSystem;
  int GroundZone;
  float ProjectionParameters[15];
  int PlaneUnitOfMeasure;
  int ElevationUnitOfMeasure;
  int PolygonSize;
  float GroundCoords[4][2];
  float ElevationBounds[2];
  float LocalRotation;
  int AccuracyCode;
  float SpatialResolution[3];
  int ProfileDimension[2];
  int ProfileSeekOffset;
  void ComputeExtentOriginAndSpacing (int extent[6], float origin[6], float spacing[6]);
  int ReadTypeARecord ();
  int ReadProfiles (vtkImageData *data);
  void ExecuteData(vtkDataObject *out);
private:
  vtkDEMReader(const vtkDEMReader&);  // Not implemented.
  void operator=(const vtkDEMReader&);  // Not implemented.
};

#endif

