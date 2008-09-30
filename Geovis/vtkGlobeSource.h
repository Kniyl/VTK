/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkGlobeSource.h

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

// .NAME vtkGlobeSource - Sphere patch with Lat/Long scalar array.
// .SECTION Description
// vtkGlobeSource will generate any "rectangular" patch of the globe
// given its Longitude-Latitude extent.  It adds two point scalar arrays
// Longitude and Latitude to the output.  These arrays can be transformed
// to generate texture coordinates for any texture map.  This source is 
// imperfect near the poles as implmented.  It should really reduce the
// longitude resolution as the triangles become slivers.

#ifndef __vtkGlobeSource_h
#define __vtkGlobeSource_h

#include "vtkPolyDataAlgorithm.h"

class vtkCellArray;
class vtkFloatArray;


class VTK_GEOVIS_EXPORT vtkGlobeSource : public vtkPolyDataAlgorithm 
{
public:
  vtkTypeRevisionMacro(vtkGlobeSource,vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);


  // Descrpition:
  // This world point will be shifted to 0,0,0.
  // Used to avoid picking bug caused by rendering errors with large offsets.
  vtkSetVector3Macro(Origin, double);

  // Description:
  // Longitude Latitude clamps.
  vtkSetClampMacro(StartLongitude,double,-180.0,180.0);
  vtkSetClampMacro(EndLongitude,double,-180.0,180.0);
  vtkSetClampMacro(StartLatitude,double,-90.0,90.0);
  vtkSetClampMacro(EndLatitude,double,-90.0,90.0);

  // Description:
  // Set the number of points in the longitude direction (ranging from
  // StartLongitude to EndLongitude).
  vtkSetClampMacro(LongitudeResolution,int,3,100);
  vtkGetMacro(LongitudeResolution,int);

  // Description:
  // Set the number of points in the latitude direction (ranging
  // from StartLatitude to EndLatitude).
  vtkSetClampMacro(LatitudeResolution,int,3,100);
  vtkGetMacro(LatitudeResolution,int);

  // Description:
  // Set radius of sphere. Default is 6356750.0
  vtkSetClampMacro(Radius,double,0.0,VTK_DOUBLE_MAX);
  vtkGetMacro(Radius,double);

  vtkSetClampMacro(CurtainHeight,double,0.0,VTK_DOUBLE_MAX);
  vtkGetMacro(CurtainHeight,double);

  // Description:
  // Cause the sphere to be tessellated with edges along the latitude
  // and longitude lines. If off, triangles are generated at non-polar
  // regions, which results in edges that are not parallel to latitude and
  // longitude lines. If on, quadrilaterals are generated everywhere
  // except at the poles. This can be useful for generating a wireframe
  // sphere with natural latitude and longitude lines.
  vtkSetMacro(QuadrilateralTessellation,int);
  vtkGetMacro(QuadrilateralTessellation,int);
  vtkBooleanMacro(QuadrilateralTessellation,int);

  // Description:
  // Construct sphere with radius=0.5 and default resolution 8 in both latitude
  // and longitude directions. longitude ranges from (-180,180) 
  // and latitude (-90,90) degrees.
  static vtkGlobeSource *New();

  // Description:
  // Calculates the normal and point on a sphere with a specified radius
  // at the spherical coordinates theta and phi.
  static void ComputeGlobePoint(
    double theta, double phi, double radius, double* point, double* normal = 0);

  // Description:
  // Calculates the spherical coordinates theta and phi based on the
  // point on a sphere.
  static void ComputeLatitudeLongitude(
    double* x, double& theta, double& phi);
  
protected:
  vtkGlobeSource();
  ~vtkGlobeSource() {}

  int RequestData(
    vtkInformation *, 
    vtkInformationVector **, 
    vtkInformationVector *);
  int RequestInformation(
    vtkInformation *, 
    vtkInformationVector **, 
    vtkInformationVector *);

  void AddPoint(
    double theta, double phi, double radius,
    vtkPoints* newPoints, vtkFloatArray* newNormals,
    vtkFloatArray* newLongitudeArray, vtkFloatArray* newLatitudeArray);


  double Origin[3];
  double Radius;
  double CurtainHeight;
  int LongitudeResolution;
  int LatitudeResolution;
  double StartLongitude;
  double EndLongitude;
  double StartLatitude;
  double EndLatitude;
  int QuadrilateralTessellation;

private:
  vtkGlobeSource(const vtkGlobeSource&);  // Not implemented.
  void operator=(const vtkGlobeSource&);  // Not implemented.
};

#endif
