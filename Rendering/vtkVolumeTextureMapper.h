/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkVolumeTextureMapper.h
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
// .NAME vtkVolumeTextureMapper - Abstract class for a volume mapper

// .SECTION Description
// vtkVolumeTextureMapper is the abstract definition of a volume mapper
// that uses a texture mapping approach.

// .SECTION see also
// vtkVolumeMapper

#ifndef __vtkVolumeTextureMapper_h
#define __vtkVolumeTextureMapper_h

#include "vtkVolumeMapper.h"
#include "vtkEncodedGradientShader.h"
#include "vtkEncodedGradientEstimator.h"

class vtkVolume;
class vtkRenderer;
class vtkRenderWindow;

class VTK_RENDERING_EXPORT vtkVolumeTextureMapper : public vtkVolumeMapper
{
public:
  vtkTypeRevisionMacro(vtkVolumeTextureMapper,vtkVolumeMapper);
  void PrintSelf( ostream& os, vtkIndent index );

  // Description:
  // Update the volume rendering pipeline by updating the scalar input
  virtual void Update();

  // Description:
  // Set / Get the gradient estimator used to estimate normals
  void SetGradientEstimator( vtkEncodedGradientEstimator *gradest );
  vtkGetObjectMacro( GradientEstimator, vtkEncodedGradientEstimator );

  // Description:
  // Get the gradient shader.
  vtkGetObjectMacro( GradientShader, vtkEncodedGradientShader );

//BTX
  // Description:
  // Allow access to the arrays / variables from the templated functions in the
  // subclasses.
  float *GetGradientOpacityArray(){return this->GradientOpacityArray;};
  unsigned char *GetRGBAArray(){return this->RGBAArray;};
  float *GetRedDiffuseShadingTable(){return this->RedDiffuseShadingTable;};
  float *GetGreenDiffuseShadingTable(){return this->GreenDiffuseShadingTable;};
  float *GetBlueDiffuseShadingTable(){return this->BlueDiffuseShadingTable;};
  float *GetRedSpecularShadingTable(){return this->RedSpecularShadingTable;};
  float *GetGreenSpecularShadingTable(){return this->GreenSpecularShadingTable;};
  float *GetBlueSpecularShadingTable(){return this->BlueSpecularShadingTable;};
  unsigned short *GetEncodedNormals(){return this->EncodedNormals;};
  unsigned char *GetGradientMagnitudes(){return this->GradientMagnitudes;};
  vtkGetMacro( Shade, int );
  vtkGetObjectMacro( RenderWindow, vtkRenderWindow );
  vtkGetVectorMacro( DataOrigin, float, 3 );
  vtkGetVectorMacro( DataSpacing, float, 3 );

  // Description:
  // WARNING: INTERNAL METHOD - NOT INTENDED FOR GENERAL USE
  // DO NOT USE THIS METHOD OUTSIDE OF THE RENDERING PROCESS
  // Render the volume
  virtual void Render(vtkRenderer *ren, vtkVolume *vol)=0;

  // Description:
  // WARNING: INTERNAL METHOD - NOT INTENDED FOR GENERAL USE
  // Values needed by the volume
  virtual float GetGradientMagnitudeScale();
  virtual float GetGradientMagnitudeBias();
  
//ETX



protected:
  vtkVolumeTextureMapper();
  ~vtkVolumeTextureMapper();

  void InitializeRender( vtkRenderer *ren, vtkVolume *vol );

  // Objects / variables  needed for shading / gradient magnitude opacity
  vtkEncodedGradientEstimator  *GradientEstimator;
  vtkEncodedGradientShader     *GradientShader;
  int                          Shade;

  float          *GradientOpacityArray;
  unsigned char  *RGBAArray;
  int            ArraySize;

  float          *RedDiffuseShadingTable;
  float          *GreenDiffuseShadingTable;
  float          *BlueDiffuseShadingTable;
  float          *RedSpecularShadingTable;
  float          *GreenSpecularShadingTable;
  float          *BlueSpecularShadingTable;

  float          DataOrigin[3];
  float          DataSpacing[3];

  unsigned short *EncodedNormals;
  unsigned char  *GradientMagnitudes;

  float          SampleDistance;
  
  vtkRenderWindow *RenderWindow;
private:
  vtkVolumeTextureMapper(const vtkVolumeTextureMapper&);  // Not implemented.
  void operator=(const vtkVolumeTextureMapper&);  // Not implemented.
};


#endif


