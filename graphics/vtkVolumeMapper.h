/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkVolumeMapper.h
  Language:  C++
  Date:      $Date$
  Version:   $Revision$


Copyright (c) 1993-1998 Ken Martin, Will Schroeder, Bill Lorensen.

This software is copyrighted by Ken Martin, Will Schroeder and Bill Lorensen.
The following terms apply to all files associated with the software unless
explicitly disclaimed in individual files. This copyright specifically does
not apply to the related textbook "The Visualization Toolkit" ISBN
013199837-4 published by Prentice Hall which is covered by its own copyright.

The authors hereby grant permission to use, copy, and distribute this
software and its documentation for any purpose, provided that existing
copyright notices are retained in all copies and that this notice is included
verbatim in any distributions. Additionally, the authors grant permission to
modify this software and its documentation for any purpose, provided that
such modifications are not distributed without the explicit consent of the
authors and that existing copyright notices are retained in all copies. Some
of the algorithms implemented by this software are patented, observe all
applicable patent law.

IN NO EVENT SHALL THE AUTHORS OR DISTRIBUTORS BE LIABLE TO ANY PARTY FOR
DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES ARISING OUT
OF THE USE OF THIS SOFTWARE, ITS DOCUMENTATION, OR ANY DERIVATIVES THEREOF,
EVEN IF THE AUTHORS HAVE BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

THE AUTHORS AND DISTRIBUTORS SPECIFICALLY DISCLAIM ANY WARRANTIES, INCLUDING,
BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
PARTICULAR PURPOSE, AND NON-INFRINGEMENT.  THIS SOFTWARE IS PROVIDED ON AN
"AS IS" BASIS, AND THE AUTHORS AND DISTRIBUTORS HAVE NO OBLIGATION TO PROVIDE
MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.


=========================================================================*/
// .NAME vtkVolumeMapper - Abstract class for a volume mapper

// .SECTION Description
// vtkVolumeMapper is the abstract definition of a volume mapper.
// Several basic types of volume mappers are supported. There are
// ray casters, which produce an image that must be merged with
// geometry, there are hardware methods that blend with geometry,
// and some combinations of these.

// .SECTION see also
// vtkVolumeRayCastMapper

#ifndef __vtkVolumeMapper_h
#define __vtkVolumeMapper_h

#include "vtkObject.h"
#include "vtkStructuredPoints.h"
#include "vtkImageToStructuredPoints.h"
#include "vtkImageCache.h"

class vtkRenderer;
class vtkVolume;

class VTK_EXPORT vtkVolumeMapper : public vtkObject
{
public:

// Description:
// Construct a vtkVolumeMapper with empty scalar input and clipping off.
  vtkVolumeMapper();

  const char *GetClassName() {return "vtkVolumeMapper";};

// Description:
// Print the vtkVolumeMapper
  void PrintSelf( ostream& os, vtkIndent index );


  virtual void Render(vtkRenderer *ren, vtkVolume *vol) = 0;

  virtual void Update();

  // Description:
  // Get the bounds of the scalar data.
  virtual float *GetBounds();

  // Description:
  // Will the hardware color and z buffers be destroyed during a render?
  virtual int DestroyHardwareBuffer( void ) = 0;

  // Description:
  // Will the image be in hardware when the render is complete?
  virtual int ImageLocatedInHardware( void ) = 0;

  // Description:
  // Get the z buffer data for the image.
  virtual float *GetZbufferData( void ) = 0;

  // Description:
  // Get the RGBA color buffer data for the image.
  virtual float *GetRGBAPixelData( void ) = 0;

  // Description:
  // Turn On/Off orthogonal clipping. (Clipping planes are
  // perpendicular to the coordinate axes.)
  vtkSetMacro(Clipping,int);
  vtkGetMacro(Clipping,int);
  vtkBooleanMacro(Clipping,int);

  float GetXminClipPlane( void ) { return this->ClippingPlanes[0]; };
  float GetXmaxClipPlane( void ) { return this->ClippingPlanes[1]; };
  float GetYminClipPlane( void ) { return this->ClippingPlanes[2]; };
  float GetYmaxClipPlane( void ) { return this->ClippingPlanes[3]; };
  float GetZminClipPlane( void ) { return this->ClippingPlanes[4]; };
  float GetZmaxClipPlane( void ) { return this->ClippingPlanes[5]; };

  // Description:
  // Set/Get the ClippingPlanes ( xmin, xmax, ymin, ymax, zmin, zmax )
  void SetClippingPlanes( float a, float b, float c, 
                          float d, float e, float f );
  void SetClippingPlanes( float p[6] ); 
  float *GetClippingPlanes( void ) { return this->ClippingPlanes; };

  // Description:
  // Set/Get the scalar input data
  vtkSetObjectMacro( ScalarInput, vtkStructuredPoints );
  void SetScalarInput(vtkImageCache *cache)
    {this->SetScalarInput(cache->GetImageToStructuredPoints()->GetOutput());}
  virtual vtkStructuredPoints *GetScalarInput() {return this->ScalarInput;};

protected:
  vtkStructuredPoints  *ScalarInput;
  int                  Clipping;
  float                ClippingPlanes[6];
  vtkTimeStamp         BuildTime;
};

inline void vtkVolumeMapper::SetClippingPlanes( 
                     float a, float b, float c, float d, float e, float f )
{
  this->ClippingPlanes[0] = a;
  this->ClippingPlanes[1] = b;
  this->ClippingPlanes[2] = c;
  this->ClippingPlanes[3] = d;
  this->ClippingPlanes[4] = e;
  this->ClippingPlanes[5] = f;
}

inline void vtkVolumeMapper::SetClippingPlanes( float p[6] )
{
  this->ClippingPlanes[0] = p[0];
  this->ClippingPlanes[1] = p[1];
  this->ClippingPlanes[2] = p[2];
  this->ClippingPlanes[3] = p[3];
  this->ClippingPlanes[4] = p[4];
  this->ClippingPlanes[5] = p[5];
}

#endif


