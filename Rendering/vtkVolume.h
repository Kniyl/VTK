/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkVolume.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkVolume - represents a volume (data & properties) in a rendered scene
//
// .SECTION Description
// vtkVolume is used to represent a volumetric entity in a rendering scene.
// It inherits functions related to the volume's position, orientation and
// origin from vtkProp3D. The volume maintains a reference to the
// volumetric data (i.e., the volume mapper). The volume also contains a
// reference to a volume property which contains all common volume rendering 
// parameters.

// .SECTION see also
// vtkAbstractVolumeMapper vtkVolumeProperty vtkProp3D

#ifndef __vtkVolume_h
#define __vtkVolume_h

#include "vtkProp3D.h"

class vtkRenderer;
class vtkPropCollection;
class vtkVolumeCollection;
class vtkWindow;
class vtkVolumeProperty;
class vtkAbstractVolumeMapper;

class VTK_RENDERING_EXPORT vtkVolume : public vtkProp3D
{
public:
  vtkTypeRevisionMacro(vtkVolume,vtkProp3D);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Creates a Volume with the following defaults: origin(0,0,0) 
  // position=(0,0,0) scale=1 visibility=1 pickable=1 dragable=1
  // orientation=(0,0,0).
  static vtkVolume *New();

  // Description:
  // Set/Get the volume mapper.
  void SetMapper(vtkAbstractVolumeMapper *mapper);
  vtkGetObjectMacro(Mapper, vtkAbstractVolumeMapper);

  // Description:
  // Set/Get the volume property.
  void SetProperty(vtkVolumeProperty *property);
  vtkVolumeProperty *GetProperty();

  // Description: 
  // For some exporters and other other operations we must be
  // able to collect all the actors or volumes. This method
  // is used in that process.
  void GetVolumes(vtkPropCollection *vc);

  // Description:
  // Update the volume rendering pipeline by updating the volume mapper
  void Update();

  // Description:
  // Get the bounds - either all six at once 
  // (xmin, xmax, ymin, ymax, zmin, zmax) or one at a time.
  float *GetBounds();
  void GetBounds(float bounds[6]) { this->vtkProp3D::GetBounds( bounds ); };
  float GetMinXBound();
  float GetMaxXBound();
  float GetMinYBound();
  float GetMaxYBound();
  float GetMinZBound();
  float GetMaxZBound();

  // Description:
  // Return the MTime also considering the property etc.
  unsigned long int GetMTime();

  // Description:
  // Return the mtime of anything that would cause the rendered image to 
  // appear differently. Usually this involves checking the mtime of the 
  // prop plus anything else it depends on such as properties, mappers,
  // etc.
  unsigned long GetRedrawMTime();

  // Description:
  // Shallow copy of this vtkVolume. Overloads the virtual vtkProp method.
  void ShallowCopy(vtkProp *prop);

//BTX
  // Description:
  // WARNING: INTERNAL METHOD - NOT INTENDED FOR GENERAL USE
  // DO NOT USE THIS METHOD OUTSIDE OF THE RENDERING PROCESS
  // Support the standard render methods.
  // Depending on the mapper type, the volume may be rendered using
  // this method (FRAMEBUFFER volume such as texture mapping will
  // be rendered this way)
  int RenderTranslucentGeometry(vtkViewport *viewport);

  // Description:
  // WARNING: INTERNAL METHOD - NOT INTENDED FOR GENERAL USE
  // Release any graphics resources that are being consumed by this volume.
  // The parameter window could be used to determine which graphic
  // resources to release.
  void ReleaseGraphicsResources(vtkWindow *);

  // Description:
  // WARNING: INTERNAL METHOD - NOT INTENDED FOR GENERAL USE
  // DO NOT USE THIS METHOD OUTSIDE OF THE RENDERING PROCESS
  float *GetCorrectedScalarOpacityArray(int);
  float *GetCorrectedScalarOpacityArray()
    {return this->GetCorrectedScalarOpacityArray(0);};

  // Description:
  // WARNING: INTERNAL METHOD - NOT INTENDED FOR GENERAL USE
  // DO NOT USE THIS METHOD OUTSIDE OF THE RENDERING PROCESS
  float *GetScalarOpacityArray(int);
  float *GetScalarOpacityArray(){return this->GetScalarOpacityArray(0);};

  // Description:
  // WARNING: INTERNAL METHOD - NOT INTENDED FOR GENERAL USE
  // DO NOT USE THIS METHOD OUTSIDE OF THE RENDERING PROCESS
  float *GetGradientOpacityArray(int);
  float *GetGradientOpacityArray(){return this->GetGradientOpacityArray(0);};

  // Description:
  // WARNING: INTERNAL METHOD - NOT INTENDED FOR GENERAL USE
  // DO NOT USE THIS METHOD OUTSIDE OF THE RENDERING PROCESS
  float *GetGrayArray(int);
  float *GetGrayArray(){return this->GetGrayArray(0);};

  // Description:
  // WARNING: INTERNAL METHOD - NOT INTENDED FOR GENERAL USE
  // DO NOT USE THIS METHOD OUTSIDE OF THE RENDERING PROCESS
  float *GetRGBArray(int);
  float *GetRGBArray(){return this->GetRGBArray(0);};

  // Description:
  // WARNING: INTERNAL METHOD - NOT INTENDED FOR GENERAL USE
  // DO NOT USE THIS METHOD OUTSIDE OF THE RENDERING PROCESS
  float  GetGradientOpacityConstant(int);
  float  GetGradientOpacityConstant()
    {return this->GetGradientOpacityConstant(0);};

  // Description:
  // WARNING: INTERNAL METHOD - NOT INTENDED FOR GENERAL USE
  // DO NOT USE THIS METHOD OUTSIDE OF THE RENDERING PROCESS
  float  GetArraySize () { return static_cast<float>(this->ArraySize); };

  // Description:
  // WARNING: INTERNAL METHOD - NOT INTENDED FOR GENERAL USE
  // DO NOT USE THIS METHOD OUTSIDE OF THE RENDERING PROCESS
  void UpdateTransferFunctions( vtkRenderer *ren );

  // Description:
  // WARNING: INTERNAL METHOD - NOT INTENDED FOR GENERAL USE
  // DO NOT USE THIS METHOD OUTSIDE OF THE RENDERING PROCESS
  void UpdateScalarOpacityforSampleSize( vtkRenderer *ren, 
                                         float sample_distance );

//ETX

protected:
  vtkVolume();
  ~vtkVolume();

  vtkAbstractVolumeMapper      *Mapper;
  vtkVolumeProperty            *Property;

  // The rgb transfer function array - for unsigned char data this
  // is 256 elements, for short or unsigned short it is 65536 elements
  // This is a sample at each scalar value of the rgb transfer
  // function.  A time stamp is kept to know when it needs rebuilding
  float                *RGBArray[VTK_MAX_VRCOMP];
  vtkTimeStamp          RGBArrayMTime[VTK_MAX_VRCOMP];

  // The gray transfer function array - for unsigned char data this
  // is 256 elements, for short or unsigned short it is 65536 elements
  // This is a sample at each scalar value of the gray transfer
  // function.  A time stamp is kept to know when it needs rebuilding
  float                *GrayArray[VTK_MAX_VRCOMP];
  vtkTimeStamp          GrayArrayMTime[VTK_MAX_VRCOMP];

  // The scalar opacity transfer function array - for unsigned char data this
  // is 256 elements, for short or unsigned short it is 65536 elements
  // This is a sample at each scalar value of the opacity transfer
  // function.  A time stamp is kept to know when it needs rebuilding
  float                *ScalarOpacityArray[VTK_MAX_VRCOMP];
  vtkTimeStamp          ScalarOpacityArrayMTime[VTK_MAX_VRCOMP];

  // The corrected scalar opacity transfer function array - this is identical
  // to the opacity transfer function array when the step size is 1.
  // In other cases, it is corrected to reflect the new material thickness
  // modeled by a step size different than 1.
  float                *CorrectedScalarOpacityArray[VTK_MAX_VRCOMP];
  vtkTimeStamp          CorrectedScalarOpacityArrayMTime[VTK_MAX_VRCOMP];

  // CorrectedStepSize is the step size currently modeled by
  // CorrectedArray.  It is used to determine when the 
  // CorrectedArray needs to be updated to match SampleDistance
  // in the volume mapper.
  float                CorrectedStepSize;

  // Number of elements in the rgb, gray, and opacity transfer function arrays
  int                  ArraySize;

  // The magnitude of gradient opacity transfer function array
  float                GradientOpacityArray[VTK_MAX_VRCOMP][256];
  float                GradientOpacityConstant[VTK_MAX_VRCOMP];
  vtkTimeStamp         GradientOpacityArrayMTime[VTK_MAX_VRCOMP];

  // Function to compute screen coverage of this volume
  float ComputeScreenCoverage( vtkViewport *vp );
  
private:
  vtkVolume(const vtkVolume&);  // Not implemented.
  void operator=(const vtkVolume&);  // Not implemented.
};

#endif

