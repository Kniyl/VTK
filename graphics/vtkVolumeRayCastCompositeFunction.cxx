/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkVolumeRayCastCompositeFunction.cxx
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
#include <math.h>

#include "vtkVolumeRayCastCompositeFunction.h"
#include "vtkVolumeProperty.h"
#include "vtkVolumeRayCastMapper.h"
#include "vtkVolume.h"

#define VTK_REMAINING_OPACITY		0.02

// This is the templated function that actually casts a ray and computes
// The composite value. This version uses nearest neighbor interpolation
// and does not perform shading.
template <class T>
static void CastRay_NN_Unshaded( vtkVolumeRayCastCompositeFunction *cast_function, 
				 T *data_ptr,
				 struct VolumeRayCastRayInfoStruct *rayInfo,
				 struct VolumeRayCastVolumeInfoStruct *volumeInfo )
{
  int             value;
  unsigned char   *grad_mag_ptr = NULL;
  float           accum_red_intensity;
  float           accum_green_intensity;
  float           accum_blue_intensity;
  float           accum_intensity;
  float           remaining_opacity;
  float           opacity;
  float           gradient_opacity;
  int             loop;
  int             xinc, yinc, zinc;
  int             voxel[3];
  float           ray_position[3];
  int             prev_voxel[3];
  float           *SOTF;
  float           *CTF;
  float           *GTF;
  float           *GOTF;
  int             offset;
  int             steps_this_ray = 0;
  int             grad_op_is_constant;
  float           gradient_opacity_constant;
  int             num_steps;
  float           *ray_start, *ray_increment;

  num_steps = rayInfo->VolumeRayNumberOfSamples;
  ray_start = rayInfo->VolumeRayStart;
  ray_increment = rayInfo->VolumeRayIncrement;
 
  SOTF =  volumeInfo->Volume->GetCorrectedScalarOpacityArray();
  CTF  =  volumeInfo->Volume->GetRGBArray();
  GTF  =  volumeInfo->Volume->GetGrayArray();
  GOTF =  volumeInfo->Volume->GetGradientOpacityArray();

  // Get the gradient opacity constant. If this number is greater than
  // or equal to 0.0, then the gradient opacity transfer function is
  // a constant at that value, otherwise it is not a constant function
  gradient_opacity_constant = volumeInfo->Volume->GetGradientOpacityConstant();
  grad_op_is_constant = ( gradient_opacity_constant >= 0.0 );

  // Move the increments into local variables
  xinc = cast_function->DataIncrement[0];
  yinc = cast_function->DataIncrement[1];
  zinc = cast_function->DataIncrement[2];

  // Initialize the ray position and voxel location
  ray_position[0] = ray_start[0];
  ray_position[1] = ray_start[1];
  ray_position[2] = ray_start[2];
  voxel[0] = vtkRoundFuncMacro( ray_position[0] );
  voxel[1] = vtkRoundFuncMacro( ray_position[1] );
  voxel[2] = vtkRoundFuncMacro( ray_position[2] );

  // So far we haven't accumulated anything
  accum_intensity         = 0.0;
  accum_red_intensity     = 0.0;
  accum_green_intensity   = 0.0;
  accum_blue_intensity    = 0.0;
  remaining_opacity       = 1.0;

  // Get a pointer to the gradient magnitudes for this volume
  if ( !grad_op_is_constant )
    {
    grad_mag_ptr = cast_function->GradientMagnitudes;
    }

  // Set up the data values for the first pass through the loop
  offset = voxel[2] * zinc + voxel[1] * yinc + voxel[0];
  value = *(data_ptr + offset);
  opacity = SOTF[value];
  if ( grad_op_is_constant )
    {
    gradient_opacity = gradient_opacity_constant;
    }
  else 
    {
    gradient_opacity = GOTF[*(grad_mag_ptr + offset)];
    }
  
  // Keep track of previous voxel to know when we step into a new one
  prev_voxel[0] = voxel[0];
  prev_voxel[1] = voxel[1];
  prev_voxel[2] = voxel[2];
  
  // Two cases - we are working with a gray or RGB transfer
  // function - break them up to make it more efficient
  if ( volumeInfo->ColorChannels == 1 ) 
    {
    // For each step along the ray
    for ( loop = 0; 
	  loop < num_steps && remaining_opacity > VTK_REMAINING_OPACITY; 
	  loop++ )
      {	    
      // We've taken another step
      steps_this_ray++;
      
      // Access the value at this voxel location
      if ( prev_voxel[0] != voxel[0] ||
	   prev_voxel[1] != voxel[1] ||
	   prev_voxel[2] != voxel[2] )
	{
	offset = voxel[2] * zinc + voxel[1] * yinc + voxel[0];
	value = *(data_ptr + offset);
	opacity = SOTF[value];

	if ( opacity )
	  {
	  if ( grad_op_is_constant )
	    {
	    gradient_opacity = gradient_opacity_constant;
	    }
	  else 
	    {
	    gradient_opacity = GOTF[*(grad_mag_ptr + offset)];
	    }
	  opacity *= gradient_opacity;
	  }

	prev_voxel[0] = voxel[0];
	prev_voxel[1] = voxel[1];
	prev_voxel[2] = voxel[2];
	}
        
      // Accumulate some light intensity and opacity
      accum_red_intensity   += ( opacity * remaining_opacity * 
				 GTF[(value)] );
      remaining_opacity *= (1.0 - opacity);
      
      // Increment our position and compute our voxel location
      ray_position[0] += ray_increment[0];
      ray_position[1] += ray_increment[1];
      ray_position[2] += ray_increment[2];
      voxel[0] = vtkRoundFuncMacro( ray_position[0] );
      voxel[1] = vtkRoundFuncMacro( ray_position[1] );
      voxel[2] = vtkRoundFuncMacro( ray_position[2] );
      }
    accum_green_intensity = accum_red_intensity;
    accum_blue_intensity = accum_red_intensity;
    }
  else if ( volumeInfo->ColorChannels == 3 )
    {
    // For each step along the ray
    for ( loop = 0; 
	  loop < num_steps && remaining_opacity > VTK_REMAINING_OPACITY; 
	  loop++ )
      {	    
      // We've taken another step
      steps_this_ray++;
      
      // Access the value at this voxel location
      if ( prev_voxel[0] != voxel[0] ||
	   prev_voxel[1] != voxel[1] ||
	   prev_voxel[2] != voxel[2] )
	{
	offset = voxel[2] * zinc + voxel[1] * yinc + voxel[0];
	value = *(data_ptr + offset);
	opacity = SOTF[value];

	if ( opacity )
	  {
	  if ( grad_op_is_constant )
	    {
	    gradient_opacity = gradient_opacity_constant;
	    }
	  else 
	    {
	    gradient_opacity = GOTF[*(grad_mag_ptr + offset)];	
	    }
	  opacity *= gradient_opacity;
	  }

	prev_voxel[0] = voxel[0];
	prev_voxel[1] = voxel[1];
	prev_voxel[2] = voxel[2];
	}


      // Accumulate some light intensity and opacity
      accum_red_intensity   += ( opacity * remaining_opacity * 
				 CTF[(value)*3] );
      accum_green_intensity += ( opacity * remaining_opacity * 
				 CTF[(value)*3 + 1] );
      accum_blue_intensity  += ( opacity * remaining_opacity * 
				 CTF[(value)*3 + 2] );
      remaining_opacity *= (1.0 - opacity);
      
      // Increment our position and compute our voxel location
      ray_position[0] += ray_increment[0];
      ray_position[1] += ray_increment[1];
      ray_position[2] += ray_increment[2];
      voxel[0] = vtkRoundFuncMacro( ray_position[0] );
      voxel[1] = vtkRoundFuncMacro( ray_position[1] );
      voxel[2] = vtkRoundFuncMacro( ray_position[2] );
      }
    }

  // Cap the intensity value at 1.0
  if ( accum_red_intensity > 1.0 )
    {
    accum_red_intensity = 1.0;
    }
  if ( accum_green_intensity > 1.0 )
    {
    accum_green_intensity = 1.0;
    }
  if ( accum_blue_intensity > 1.0 )
    {
    accum_blue_intensity = 1.0;
    }
  
  if( remaining_opacity < VTK_REMAINING_OPACITY )
    {
    remaining_opacity = 0.0;
    }

  // Set the return pixel value.  The depth value is currently useless and
  // should be fixed.  What should depth be in this case?  First 
  // non-opaque or total opacity or what??
  rayInfo->RayColor[0] = accum_red_intensity;
  rayInfo->RayColor[1] = accum_green_intensity;
  rayInfo->RayColor[2] = accum_blue_intensity;
  rayInfo->RayColor[3] = 1.0 - remaining_opacity;
  rayInfo->VolumeRayStepsTaken = steps_this_ray;

  if ( remaining_opacity < 1.0 )
    {
    rayInfo->RayDepth = 0.0;
    }
  else 
    {
    rayInfo->RayDepth = VTK_LARGE_FLOAT;
    }
  
}


// This is the templated function that actually casts a ray and computes
// the composite value. This version uses nearest neighbor and does
// perform shading.
template <class T>
static void CastRay_NN_Shaded( vtkVolumeRayCastCompositeFunction *cast_function, 
			       T *data_ptr,
			       struct VolumeRayCastRayInfoStruct *rayInfo,
			       struct VolumeRayCastVolumeInfoStruct *volumeInfo )
{
  int             value = 0;
  unsigned char   *grad_mag_ptr = NULL;
  float           accum_red_intensity;
  float           accum_green_intensity;
  float           accum_blue_intensity;
  float           remaining_opacity;
  float           opacity = 0.0;
  float           gradient_opacity;
  int             loop;
  int             xinc, yinc, zinc;
  int             voxel[3];
  float           ray_position[3];
  int             prev_voxel[3];
  float           *SOTF;
  float           *CTF;
  float           *GTF;
  float           *GOTF;
  float           *red_d_shade, *green_d_shade, *blue_d_shade;
  float           *red_s_shade, *green_s_shade, *blue_s_shade;
  unsigned short  *encoded_normals;
  float           red_shaded_value   = 0.0;
  float           green_shaded_value = 0.0;
  float           blue_shaded_value  = 0.0;
  int             offset = 0;
  int             steps_this_ray = 0;
  int             grad_op_is_constant;
  float           gradient_opacity_constant;
  int             num_steps;
  float           *ray_start, *ray_increment;

  num_steps = rayInfo->VolumeRayNumberOfSamples;
  ray_start = rayInfo->VolumeRayStart;
  ray_increment = rayInfo->VolumeRayIncrement;
 
  // Get diffuse shading table pointers
  red_d_shade = volumeInfo->RedDiffuseShadingTable;
  green_d_shade = volumeInfo->GreenDiffuseShadingTable;
  blue_d_shade = volumeInfo->BlueDiffuseShadingTable;

  // Get specular shading table pointers
  red_s_shade = volumeInfo->RedSpecularShadingTable;
  green_s_shade = volumeInfo->GreenSpecularShadingTable;
  blue_s_shade = volumeInfo->BlueSpecularShadingTable;

  // Get a pointer to the encoded normals for this volume
  encoded_normals = cast_function->EncodedNormals;

  // Get the scalar opacity transfer function for this volume (which maps
  // scalar input values to opacities)
  SOTF =  volumeInfo->Volume->GetCorrectedScalarOpacityArray();

  // Get the color transfer function for this volume (which maps
  // scalar input values to RGB values)
  CTF =  volumeInfo->Volume->GetRGBArray();
  GTF =  volumeInfo->Volume->GetGrayArray();

  // Get the gradient opacity transfer function for this volume (which maps
  // gradient magnitudes to opacities)
  GOTF =  volumeInfo->Volume->GetGradientOpacityArray();

  // Get the gradient opacity constant. If this number is greater than
  // or equal to 0.0, then the gradient opacity transfer function is
  // a constant at that value, otherwise it is not a constant function
  gradient_opacity_constant = volumeInfo->Volume->GetGradientOpacityConstant();
  grad_op_is_constant = ( gradient_opacity_constant >= 0.0 );

  // Get a pointer to the gradient magnitudes for this volume
  if ( !grad_op_is_constant )
    {
    grad_mag_ptr = cast_function->GradientMagnitudes;
    }

  // Move the increments into local variables
  xinc = cast_function->DataIncrement[0];
  yinc = cast_function->DataIncrement[1];
  zinc = cast_function->DataIncrement[2];

  // Initialize the ray position and voxel location
  ray_position[0] = ray_start[0];
  ray_position[1] = ray_start[1];
  ray_position[2] = ray_start[2];

  voxel[0] = vtkRoundFuncMacro( ray_position[0] );
  voxel[1] = vtkRoundFuncMacro( ray_position[1] );
  voxel[2] = vtkRoundFuncMacro( ray_position[2] );

  // So far we haven't accumulated anything
  accum_red_intensity     = 0.0;
  accum_green_intensity   = 0.0;
  accum_blue_intensity    = 0.0;
  remaining_opacity       = 1.0;

  // Keep track of previous voxel to know when we step into a new one  
  prev_voxel[0] = voxel[0]-1;
  prev_voxel[1] = voxel[1]-1;
  prev_voxel[2] = voxel[2]-1;
  
  // Two cases - we are working with a gray or RGB transfer
  // function - break them up to make it more efficient
  if ( volumeInfo->ColorChannels == 1 ) 
    {
    // For each step along the ray
    for ( loop = 0; 
	  loop < num_steps && remaining_opacity > VTK_REMAINING_OPACITY; loop++ )
      {	    
      // We've taken another step
      steps_this_ray++;
      
      // Access the value at this voxel location and compute
      // opacity and shaded value
      if ( prev_voxel[0] != voxel[0] ||
	   prev_voxel[1] != voxel[1] ||
	   prev_voxel[2] != voxel[2] )
	{
	offset = voxel[2] * zinc + voxel[1] * yinc + voxel[0];
	value = *(data_ptr + offset);
      
	// Get the opacity contributed by the scalar opacity transfer function
	opacity = SOTF[value];

	// Multiply by the opacity contributed by the gradient magnitude
	// transfer function (don't both if opacity is already 0)
	if ( opacity )
	  {
	  if ( grad_op_is_constant )
	    {
	    gradient_opacity = gradient_opacity_constant;
	    }
	  else 
	    {
	    gradient_opacity = GOTF[*(grad_mag_ptr + offset)];
	    }
	  
	  opacity *= gradient_opacity;
	  
	  }

	// Compute the red shaded value (only if there is some opacity)
	// This is grey-scale so green and blue are the same as red
	if ( opacity )
	  {
	  red_shaded_value = opacity * remaining_opacity *
	    ( red_d_shade[*(encoded_normals + offset)] * GTF[value] +
	      red_s_shade[*(encoded_normals + offset)] );
	  }
	else
	  {
	  red_shaded_value = 0.0;
	  }

	prev_voxel[0] = voxel[0];
	prev_voxel[1] = voxel[1];
	prev_voxel[2] = voxel[2];
	}
    
    
      // Accumulate the shaded intensity and opacity of this sample
      accum_red_intensity += red_shaded_value;
      remaining_opacity *= (1.0 - opacity);
    
      // Increment our position and compute our voxel location
      ray_position[0] += ray_increment[0];
      ray_position[1] += ray_increment[1];
      ray_position[2] += ray_increment[2];
      voxel[0] = vtkRoundFuncMacro( ray_position[0] );
      voxel[1] = vtkRoundFuncMacro( ray_position[1] );
      voxel[2] = vtkRoundFuncMacro( ray_position[2] );
      }
    accum_green_intensity = accum_red_intensity;
    accum_blue_intensity = accum_red_intensity;
    }
  else if ( volumeInfo->ColorChannels == 3 )
    {
    // For each step along the ray
    for ( loop = 0; 
	  loop < num_steps && remaining_opacity > VTK_REMAINING_OPACITY; loop++ )
      {	    
      // We've taken another step
      steps_this_ray++;
      
      // Access the value at this voxel location and compute
      // opacity and shaded value
      if ( prev_voxel[0] != voxel[0] ||
	   prev_voxel[1] != voxel[1] ||
	   prev_voxel[2] != voxel[2] )
	{
	offset = voxel[2] * zinc + voxel[1] * yinc + voxel[0];
	value = *(data_ptr + offset);
      
	// Get the opacity contributed by the scalar opacity transfer function
	opacity = SOTF[value];

	// Multiply by the opacity contributed by the gradient magnitude
	// transfer function (don't both if opacity is already 0)
	if ( opacity )
	  {
	  if ( grad_op_is_constant )
	    {
	    gradient_opacity = gradient_opacity_constant;
	    }
	  else 
	    {
	    gradient_opacity = GOTF[*(grad_mag_ptr + offset)];
	    }
	  
	  opacity *= gradient_opacity;	  
	  }	

	// Compute the red, green, and blue shaded value (only if there
	// is some opacity)
	if ( opacity )
	  {
	  red_shaded_value = opacity *  remaining_opacity *
	    ( red_d_shade[*(encoded_normals + offset)] * CTF[value*3] +
	      red_s_shade[*(encoded_normals + offset)] );
	  green_shaded_value = opacity *  remaining_opacity *
	    ( green_d_shade[*(encoded_normals + offset)] * CTF[value*3 + 1] +
	      green_s_shade[*(encoded_normals + offset)] );
	  blue_shaded_value = opacity *  remaining_opacity *
	    ( blue_d_shade[*(encoded_normals + offset)] * CTF[value*3 + 2] +
	      blue_s_shade[*(encoded_normals + offset)] );
	  }
	else
	  {
	  red_shaded_value = 0.0;
	  green_shaded_value = 0.0;
	  blue_shaded_value = 0.0;
	  }

	prev_voxel[0] = voxel[0];
	prev_voxel[1] = voxel[1];
	prev_voxel[2] = voxel[2];
	}
    
    
      // Accumulate the shaded intensity and opacity of this sample
      accum_red_intensity += red_shaded_value;
      accum_green_intensity += green_shaded_value;
      accum_blue_intensity += blue_shaded_value;
      remaining_opacity *= (1.0 - opacity);
    
      // Increment our position and compute our voxel location
      ray_position[0] += ray_increment[0];
      ray_position[1] += ray_increment[1];
      ray_position[2] += ray_increment[2];
      voxel[0] = vtkRoundFuncMacro( ray_position[0] );
      voxel[1] = vtkRoundFuncMacro( ray_position[1] );
      voxel[2] = vtkRoundFuncMacro( ray_position[2] );
      }
    }
  
  // Cap the intensities at 1.0
  if ( accum_red_intensity > 1.0 )
    {
    accum_red_intensity = 1.0;
    }
  if ( accum_green_intensity > 1.0 )
    {
    accum_green_intensity = 1.0;
    }
  if ( accum_blue_intensity > 1.0 )
    {
    accum_blue_intensity = 1.0;
    }
  
  if( remaining_opacity < VTK_REMAINING_OPACITY )
    {
    remaining_opacity = 0.0;
    }
  
  // Set the return pixel value.  The depth value is currently useless and
  // should be fixed.  What should depth be in this case?  First 
  // non-opaque or total opacity or what??
  rayInfo->RayColor[0] = accum_red_intensity;
  rayInfo->RayColor[1] = accum_green_intensity;
  rayInfo->RayColor[2] = accum_blue_intensity;
  rayInfo->RayColor[3] = 1.0 - remaining_opacity;
  rayInfo->VolumeRayStepsTaken = steps_this_ray;
  
  if ( remaining_opacity < 1.0 )
    {
    rayInfo->RayDepth = 0.0;
    }
  else 
    {
    rayInfo->RayDepth = VTK_LARGE_FLOAT;
    }
}

// This is the templated function that actually casts a ray and computes
// the composite value.  This version uses trilinear interpolation and
// does not compute shading
template <class T>
static void CastRay_TrilinSample_Unshaded( 
					  vtkVolumeRayCastCompositeFunction *cast_function, 
					  T *data_ptr,
					  struct VolumeRayCastRayInfoStruct *rayInfo,
					  struct VolumeRayCastVolumeInfoStruct *volumeInfo )
{
  unsigned char   *grad_mag_ptr = NULL;
  unsigned char   *gmptr = NULL;
  float           accum_intensity;
  float           accum_red_intensity;
  float           accum_green_intensity;
  float           accum_blue_intensity;
  float           remaining_opacity;
  float           red_value, green_value, blue_value;
  float           opacity;
  int             loop;
  int             xinc, yinc, zinc;
  int             voxel[3];
  float           ray_position[3];
  float           A, B, C, D, E, F, G, H;
  int             Binc, Cinc, Dinc, Einc, Finc, Ginc, Hinc;
  T               *dptr;
  float           *SOTF;
  float           *CTF;
  float           *GTF;
  float           *GOTF;
  float           x, y, z, t1, t2, t3;
  int             offset;
  int             steps_this_ray = 0;
  float           gradient_value;
  float           scalar_value;
  int             grad_op_is_constant;
  float           gradient_opacity_constant;
  int             num_steps;
  float           *ray_start, *ray_increment;

  num_steps = rayInfo->VolumeRayNumberOfSamples;
  ray_start = rayInfo->VolumeRayStart;
  ray_increment = rayInfo->VolumeRayIncrement;

  // Get the scalar opacity transfer function which maps scalar input values
  // to opacities
  SOTF =  volumeInfo->Volume->GetCorrectedScalarOpacityArray();

  // Get the color transfer function which maps scalar input values
  // to RGB colors
  CTF =  volumeInfo->Volume->GetRGBArray();
  GTF =  volumeInfo->Volume->GetGrayArray();

  // Get the gradient opacity transfer function for this volume (which maps
  // gradient magnitudes to opacities)
  GOTF =  volumeInfo->Volume->GetGradientOpacityArray();

  // Get the gradient opacity constant. If this number is greater than
  // or equal to 0.0, then the gradient opacity transfer function is
  // a constant at that value, otherwise it is not a constant function
  gradient_opacity_constant = volumeInfo->Volume->GetGradientOpacityConstant();
  grad_op_is_constant = ( gradient_opacity_constant >= 0.0 );

  // Get a pointer to the gradient magnitudes for this volume
  if ( !grad_op_is_constant )
    {
    grad_mag_ptr = cast_function->GradientMagnitudes;
    }

  // Move the increments into local variables
  xinc = cast_function->DataIncrement[0];
  yinc = cast_function->DataIncrement[1];
  zinc = cast_function->DataIncrement[2];

  // Initialize the ray position and voxel location
  ray_position[0] = ray_start[0];
  ray_position[1] = ray_start[1];
  ray_position[2] = ray_start[2];
  voxel[0] = (int)( ray_position[0] );
  voxel[1] = (int)( ray_position[1] );
  voxel[2] = (int)( ray_position[2] );

  // So far we have not accumulated anything
  accum_intensity         = 0.0;
  accum_red_intensity     = 0.0;
  accum_green_intensity   = 0.0;
  accum_blue_intensity    = 0.0;
  remaining_opacity       = 1.0;

  // Compute the increments to get to the other 7 voxel vertices from A
  Binc = xinc;
  Cinc = yinc;
  Dinc = xinc + yinc;
  Einc = zinc;
  Finc = zinc + xinc;
  Ginc = zinc + yinc;
  Hinc = zinc + xinc + yinc;
  
  // Two cases - we are working with a gray or RGB transfer
  // function - break them up to make it more efficient
  if ( volumeInfo->ColorChannels == 1 ) 
    {
    // For each step along the ray
    for ( loop = 0; 
	  loop < num_steps && remaining_opacity > VTK_REMAINING_OPACITY; 
	  loop++ )
      {	    
      // We've taken another step
      steps_this_ray++;
      
      offset = voxel[2] * zinc + voxel[1] * yinc + voxel[0];
      dptr = data_ptr + offset;
      
      A = *(dptr);
      B = *(dptr + Binc);
      C = *(dptr + Cinc);
      D = *(dptr + Dinc);
      E = *(dptr + Einc);
      F = *(dptr + Finc);
      G = *(dptr + Ginc);
      H = *(dptr + Hinc);
      
      // Compute our offset in the voxel, and use that to trilinearly
      // interpolate the value
      x = ray_position[0] - (float) voxel[0];
      y = ray_position[1] - (float) voxel[1];
      z = ray_position[2] - (float) voxel[2];
      
      t1 = 1.0 - x;
      t2 = 1.0 - y;
      t3 = 1.0 - z;
      
      scalar_value = 
	A * t1 * t2 * t3 +
	B *  x * t2 * t3 +
	C * t1 *  y * t3 + 
	D *  x *  y * t3 +
	E * t1 * t2 *  z + 
	F *  x * t2 *  z + 
	G * t1 *  y *  z + 
	H *  x *  y *  z;
      
      if ( scalar_value < 0.0 ) 
	{
	scalar_value = 0.0;
	}
      else if ( scalar_value > volumeInfo->Volume->GetArraySize() - 1 )
	{
	scalar_value = volumeInfo->Volume->GetArraySize() - 1;
	}
      
      opacity = SOTF[(int)scalar_value];
      
      if ( opacity )
	{
	if ( !grad_op_is_constant )
	  {
	  gmptr = grad_mag_ptr + offset;
      
	  A = *(gmptr);
	  B = *(gmptr + Binc);
	  C = *(gmptr + Cinc);
	  D = *(gmptr + Dinc);
	  E = *(gmptr + Einc);
	  F = *(gmptr + Finc);
	  G = *(gmptr + Ginc);
	  H = *(gmptr + Hinc);
	  
	  gradient_value = 
	    A * t1 * t2 * t3 +
	    B *  x * t2 * t3 +
	    C * t1 *  y * t3 + 
	    D *  x *  y * t3 +
	    E * t1 * t2 *  z + 
	    F *  x * t2 *  z + 
	    G * t1 *  y *  z + 
	    H *  x *  y *  z;
      
	  if ( gradient_value < 0.0 )
	    {
	    gradient_value = 0.0;
	    }
	  else if ( gradient_value > 255.0 )
	    {
	    gradient_value = 255.0;
	    }

	  opacity *= GOTF[(int)gradient_value];
	  }
	else
	  {
	  opacity *= gradient_opacity_constant;
	  }
	red_value   = opacity * GTF[((int)scalar_value)];
	
	// Accumulate intensity and opacity for this sample location
	accum_red_intensity   += remaining_opacity * red_value;
	remaining_opacity *= (1.0 - opacity);
	}
    
      // Increment our position and compute our voxel location
      ray_position[0] += ray_increment[0];
      ray_position[1] += ray_increment[1];
      ray_position[2] += ray_increment[2];      
      voxel[0] = (int)( ray_position[0] );
      voxel[1] = (int)( ray_position[1] );
      voxel[2] = (int)( ray_position[2] );
      }
    accum_green_intensity = accum_red_intensity;
    accum_blue_intensity = accum_red_intensity;
    }
  else if ( volumeInfo->ColorChannels == 3 )
    {
    // For each step along the ray
    for ( loop = 0; 
	  loop < num_steps && remaining_opacity > VTK_REMAINING_OPACITY; 
	  loop++ )
      {	    
      // We've taken another step
      steps_this_ray++;
      
      offset = voxel[2] * zinc + voxel[1] * yinc + voxel[0];
      dptr = data_ptr + offset;
      
      A = *(dptr);
      B = *(dptr + Binc);
      C = *(dptr + Cinc);
      D = *(dptr + Dinc);
      E = *(dptr + Einc);
      F = *(dptr + Finc);
      G = *(dptr + Ginc);
      H = *(dptr + Hinc);
      
      // Compute our offset in the voxel, and use that to trilinearly
      // interpolate the value
      x = ray_position[0] - (float) voxel[0];
      y = ray_position[1] - (float) voxel[1];
      z = ray_position[2] - (float) voxel[2];
      
      t1 = 1.0 - x;
      t2 = 1.0 - y;
      t3 = 1.0 - z;
      
      scalar_value = 
	A * t1 * t2 * t3 +
	B *  x * t2 * t3 +
	C * t1 *  y * t3 + 
	D *  x *  y * t3 +
	E * t1 * t2 *  z + 
	F *  x * t2 *  z + 
	G * t1 *  y *  z + 
	H *  x *  y *  z;
      
      if ( scalar_value < 0.0 ) 
	{
	scalar_value = 0.0;
	}
      else if ( scalar_value > volumeInfo->Volume->GetArraySize() - 1 )
	{
	scalar_value = volumeInfo->Volume->GetArraySize() - 1;
	}
      
      opacity = SOTF[(int)scalar_value];
      
      if ( opacity )
	{
	if ( !grad_op_is_constant )
	  {
	  gmptr = grad_mag_ptr + offset;
	  
	  A = *(gmptr);
	  B = *(gmptr + Binc);
	  C = *(gmptr + Cinc);
	  D = *(gmptr + Dinc);
	  E = *(gmptr + Einc);
	  F = *(gmptr + Finc);
	  G = *(gmptr + Ginc);
	  H = *(gmptr + Hinc);
	  
	  gradient_value = 
	    A * t1 * t2 * t3 +
	    B *  x * t2 * t3 +
	    C * t1 *  y * t3 + 
	    D *  x *  y * t3 +
	    E * t1 * t2 *  z + 
	    F *  x * t2 *  z + 
	    G * t1 *  y *  z + 
	    H *  x *  y *  z;
	  
	  if ( gradient_value < 0.0 )
	    {
	    gradient_value = 0.0;
	    }
	  else if ( gradient_value > 255.0 )
	    {
	    gradient_value = 255.0;
	    }

	  opacity *= GOTF[(int)gradient_value];
	  }
	else
	  {
	  opacity *= gradient_opacity_constant;
	  }

	red_value   = opacity * CTF[((int)scalar_value) * 3    ];
	green_value = opacity * CTF[((int)scalar_value) * 3 + 1];
	blue_value  = opacity * CTF[((int)scalar_value) * 3 + 2];
	
	// Accumulate intensity and opacity for this sample location
	accum_red_intensity   += remaining_opacity * red_value;
	accum_green_intensity += remaining_opacity * green_value;
	accum_blue_intensity  += remaining_opacity * blue_value;
	remaining_opacity *= (1.0 - opacity);
	}
    
      // Increment our position and compute our voxel location
      ray_position[0] += ray_increment[0];
      ray_position[1] += ray_increment[1];
      ray_position[2] += ray_increment[2];      
      voxel[0] = (int)( ray_position[0] );
      voxel[1] = (int)( ray_position[1] );
      voxel[2] = (int)( ray_position[2] );
      }
    }

  // Cap the intensity value at 1.0
  if ( accum_red_intensity > 1.0 )
    {
    accum_red_intensity = 1.0;
    }
  if ( accum_green_intensity > 1.0 )
    {
    accum_green_intensity = 1.0;
    }
  if ( accum_blue_intensity > 1.0 )
    {
    accum_blue_intensity = 1.0;
    }
  
  if( remaining_opacity < VTK_REMAINING_OPACITY )
    {
    remaining_opacity = 0.0;
    }

  // Set the return pixel value.  The depth value is currently useless and
  // should be fixed.  What should depth be in this case?  First 
  // non-opaque or total opacity or what??
  rayInfo->RayColor[0] = accum_red_intensity;
  rayInfo->RayColor[1] = accum_green_intensity;
  rayInfo->RayColor[2] = accum_blue_intensity;
  rayInfo->RayColor[3] = 1.0 - remaining_opacity;
  rayInfo->VolumeRayStepsTaken = steps_this_ray;

  if ( remaining_opacity < 1.0 )
    {
    rayInfo->RayDepth = 0.0;
    }
  else 
    {
    rayInfo->RayDepth = VTK_LARGE_FLOAT;
    }
}

// This is the templated function that actually casts a ray and computes
// the composite value.  This version uses trilinear interpolation, and
// does perform shading.
template <class T>
static void CastRay_TrilinSample_Shaded( 
					vtkVolumeRayCastCompositeFunction *cast_function, 
					T *data_ptr,
					struct VolumeRayCastRayInfoStruct *rayInfo,
					struct VolumeRayCastVolumeInfoStruct *volumeInfo )
{
  unsigned char   *grad_mag_ptr = NULL;
  unsigned char   *gmptr = NULL;
  float           accum_red_intensity;
  float           accum_green_intensity;
  float           accum_blue_intensity;
  float           remaining_opacity;
  float           opacity;
  int             loop;
  int             xinc, yinc, zinc;
  int             voxel[3];
  float           ray_position[3];
  float           A, B, C, D, E, F, G, H;
  int             A_n, B_n, C_n, D_n, E_n, F_n, G_n, H_n;
  float           final_rd, final_gd, final_bd;
  float           final_rs, final_gs, final_bs;
  int             Binc, Cinc, Dinc, Einc, Finc, Ginc, Hinc;
  T               *dptr;
  float           *SOTF;
  float           *CTF;
  float           *GTF;
  float           *GOTF;
  float           x, y, z, t1, t2, t3;
  float           tA, tB, tC, tD, tE, tF, tG, tH;
  float           *red_d_shade, *green_d_shade, *blue_d_shade;
  float           *red_s_shade, *green_s_shade, *blue_s_shade;
  unsigned short  *encoded_normals, *nptr;
  float           red_shaded_value, green_shaded_value, blue_shaded_value;
  int             offset;
  int             steps_this_ray = 0;
  int             gradient_value;
  int             scalar_value;
  float           r, g, b;
  int             grad_op_is_constant;
  float           gradient_opacity_constant;
  int             num_steps;
  float           *ray_start, *ray_increment;

  num_steps = rayInfo->VolumeRayNumberOfSamples;
  ray_start = rayInfo->VolumeRayStart;
  ray_increment = rayInfo->VolumeRayIncrement;

  // Get diffuse shading table pointers
  red_d_shade = volumeInfo->RedDiffuseShadingTable;
  green_d_shade = volumeInfo->GreenDiffuseShadingTable;
  blue_d_shade = volumeInfo->BlueDiffuseShadingTable;


  // Get diffuse shading table pointers
  red_s_shade = volumeInfo->RedSpecularShadingTable;
  green_s_shade = volumeInfo->GreenSpecularShadingTable;
  blue_s_shade = volumeInfo->BlueSpecularShadingTable;

  // Get a pointer to the encoded normals for this volume
  encoded_normals = cast_function->EncodedNormals;

  // Get the scalar opacity transfer function which maps scalar input values
  // to opacities
  SOTF =  volumeInfo->Volume->GetCorrectedScalarOpacityArray();

  // Get the color transfer function which maps scalar input values
  // to RGB values
  CTF =  volumeInfo->Volume->GetRGBArray();
  GTF =  volumeInfo->Volume->GetGrayArray();

  // Get the gradient opacity transfer function for this volume (which maps
  // gradient magnitudes to opacities)
  GOTF =  volumeInfo->Volume->GetGradientOpacityArray();

  // Get the gradient opacity constant. If this number is greater than
  // or equal to 0.0, then the gradient opacity transfer function is
  // a constant at that value, otherwise it is not a constant function
  gradient_opacity_constant = volumeInfo->Volume->GetGradientOpacityConstant();
  grad_op_is_constant = ( gradient_opacity_constant >= 0.0 );

  // Get a pointer to the gradient magnitudes for this volume
  if ( !grad_op_is_constant )
    {
    grad_mag_ptr = cast_function->GradientMagnitudes;
    }

  // Move the increments into local variables
  xinc = cast_function->DataIncrement[0];
  yinc = cast_function->DataIncrement[1];
  zinc = cast_function->DataIncrement[2];

  // Initialize the ray position and voxel location
  ray_position[0] = ray_start[0];
  ray_position[1] = ray_start[1];
  ray_position[2] = ray_start[2];
  voxel[0] = (int)( ray_position[0] );
  voxel[1] = (int)( ray_position[1] );
  voxel[2] = (int)( ray_position[2] );

  // So far we haven't accumulated anything
  accum_red_intensity   = 0.0;
  accum_green_intensity = 0.0;
  accum_blue_intensity  = 0.0;
  remaining_opacity     = 1.0;

  // Compute the increments to get to the other 7 voxel vertices from A
  Binc = xinc;
  Cinc = yinc;
  Dinc = xinc + yinc;
  Einc = zinc;
  Finc = zinc + xinc;
  Ginc = zinc + yinc;
  Hinc = zinc + xinc + yinc;
  
  // Two cases - we are working with a gray or RGB transfer
  // function - break them up to make it more efficient
  if ( volumeInfo->ColorChannels == 1 ) 
    {
    // For each step along the ray
    for ( loop = 0; 
	  loop < num_steps && remaining_opacity > VTK_REMAINING_OPACITY; 
	  loop++ )
      {	    
      // We've taken another step
      steps_this_ray++;
      
      offset = voxel[2] * zinc + voxel[1] * yinc + voxel[0];
      dptr = data_ptr + offset;
      nptr = encoded_normals + offset;
    
      A = *(dptr);
      B = *(dptr + Binc);
      C = *(dptr + Cinc);
      D = *(dptr + Dinc);
      E = *(dptr + Einc);
      F = *(dptr + Finc);
      G = *(dptr + Ginc);
      H = *(dptr + Hinc);
      
      // Compute our offset in the voxel, and use that to trilinearly
      // interpolate a value
      x = ray_position[0] - (float) voxel[0];
      y = ray_position[1] - (float) voxel[1];
      z = ray_position[2] - (float) voxel[2];
      
      t1 = 1.0 - x;
      t2 = 1.0 - y;
      t3 = 1.0 - z;
      
      tA = t1 * t2 * t3;
      tB =  x * t2 * t3;
      tC = t1 *  y * t3;
      tD =  x *  y * t3;
      tE = t1 * t2 *  z;
      tF =  x * t2 *  z;
      tG = t1 *  y *  z;
      tH =  x *  y *  z;
      
      scalar_value = (int) (
			    A * tA + B * tB + C * tC + D * tD + 
			    E * tE + F * tF + G * tG + H * tH );
      
      if ( scalar_value < 0 ) 
	{
	scalar_value = 0;
	}
      else if ( scalar_value > volumeInfo->Volume->GetArraySize() - 1 )
	{
	scalar_value = volumeInfo->Volume->GetArraySize() - 1;
	}
      
      opacity = SOTF[scalar_value];

      // If we have some opacity based on the scalar value transfer function,
      // then multiply by the opacity from the gradient magnitude transfer
      // function
      if ( opacity )
	{
	if ( !grad_op_is_constant )
	  {
	  gmptr = grad_mag_ptr + offset;
      
	  A = *(gmptr);
	  B = *(gmptr + Binc);
	  C = *(gmptr + Cinc);
	  D = *(gmptr + Dinc);
	  E = *(gmptr + Einc);
	  F = *(gmptr + Finc);
	  G = *(gmptr + Ginc);
	  H = *(gmptr + Hinc);
	  
	  gradient_value = (int) (
				  A * tA + B * tB + C * tC + D * tD + 
				  E * tE + F * tF + G * tG + H * tH );
	  if ( gradient_value < 0 )
	    {
	    gradient_value = 0;
	    }
	  else if ( gradient_value > 255 )
	    {
	    gradient_value = 255;
	    }
	  
	  opacity *= GOTF[gradient_value];
	  }
	else
	  {
	  opacity *= gradient_opacity_constant;
	  }
	}

      // If we have a combined opacity value, then compute the shading
      if ( opacity )
	{
	A_n = *(nptr);
	B_n = *(nptr + Binc);
	C_n = *(nptr + Cinc);
	D_n = *(nptr + Dinc);
	E_n = *(nptr + Einc);
	F_n = *(nptr + Finc);
	G_n = *(nptr + Ginc);
	H_n = *(nptr + Hinc);

	final_rd = 
	  red_d_shade[ A_n ] * tA + red_d_shade[ B_n ] * tB + 	
	  red_d_shade[ C_n ] * tC + red_d_shade[ D_n ] * tD + 
	  red_d_shade[ E_n ] * tE + red_d_shade[ F_n ] * tF +	
	  red_d_shade[ G_n ] * tG + red_d_shade[ H_n ] * tH;
	
	final_rs = 
	  red_s_shade[ A_n ] * tA + red_s_shade[ B_n ] * tB + 	
	  red_s_shade[ C_n ] * tC + red_s_shade[ D_n ] * tD + 
	  red_s_shade[ E_n ] * tE + red_s_shade[ F_n ] * tF +	
	  red_s_shade[ G_n ] * tG + red_s_shade[ H_n ] * tH;
	
	r = GTF[(scalar_value)];

	// For this sample we have do not yet have any opacity or
	// shaded intensity yet
	red_shaded_value   = opacity * ( final_rd * r + final_rs );
	
	// Accumulate intensity and opacity for this sample location   
	accum_red_intensity   += red_shaded_value * remaining_opacity;
	remaining_opacity *= (1.0 - opacity);
	}

      // Increment our position and compute our voxel location
      ray_position[0] += ray_increment[0];
      ray_position[1] += ray_increment[1];
      ray_position[2] += ray_increment[2];      
      voxel[0] = (int)( ray_position[0] );
      voxel[1] = (int)( ray_position[1] );
      voxel[2] = (int)( ray_position[2] );
      }
    accum_green_intensity = accum_red_intensity;
    accum_blue_intensity = accum_red_intensity;
    }
  else if ( volumeInfo->ColorChannels == 3 )
    {
    // For each step along the ray
    for ( loop = 0; 
	  loop < num_steps && remaining_opacity > VTK_REMAINING_OPACITY; 
	  loop++ )
      {	    
      // We've taken another step
      steps_this_ray++;
      
      offset = voxel[2] * zinc + voxel[1] * yinc + voxel[0];
      dptr = data_ptr + offset;
      nptr = encoded_normals + offset;
    
      A = *(dptr);
      B = *(dptr + Binc);
      C = *(dptr + Cinc);
      D = *(dptr + Dinc);
      E = *(dptr + Einc);
      F = *(dptr + Finc);
      G = *(dptr + Ginc);
      H = *(dptr + Hinc);
      
      // Compute our offset in the voxel, and use that to trilinearly
      // interpolate a value
      x = ray_position[0] - (float) voxel[0];
      y = ray_position[1] - (float) voxel[1];
      z = ray_position[2] - (float) voxel[2];
      
      t1 = 1.0 - x;
      t2 = 1.0 - y;
      t3 = 1.0 - z;
      
      tA = t1 * t2 * t3;
      tB =  x * t2 * t3;
      tC = t1 *  y * t3;
      tD =  x *  y * t3;
      tE = t1 * t2 *  z;
      tF =  x * t2 *  z;
      tG = t1 *  y *  z;
      tH =  x *  y *  z;
      
      scalar_value = (int) (
			    A * tA + B * tB + C * tC + D * tD + 
			    E * tE + F * tF + G * tG + H * tH );
      
      if ( scalar_value < 0 ) 
	{
	scalar_value = 0;
	}
      else if ( scalar_value > volumeInfo->Volume->GetArraySize() - 1 )
	{
	scalar_value = volumeInfo->Volume->GetArraySize() - 1;
	}
      
      opacity = SOTF[scalar_value];
      
      if ( opacity )
	{
	  if ( !grad_op_is_constant )
	    {
	    gmptr = grad_mag_ptr + offset;
      
	    A = *(gmptr);
	    B = *(gmptr + Binc);
	    C = *(gmptr + Cinc);
	    D = *(gmptr + Dinc);
	    E = *(gmptr + Einc);
	    F = *(gmptr + Finc);
	    G = *(gmptr + Ginc);
	    H = *(gmptr + Hinc);
	    
	    gradient_value = (int) (
				    A * tA + B * tB + C * tC + D * tD + 
				    E * tE + F * tF + G * tG + H * tH );
	    if ( gradient_value < 0 )
	      {
	      gradient_value = 0;
	      }
	    else if ( gradient_value > 255 )
	      {
	      gradient_value = 255;
	      }

	    opacity *= GOTF[gradient_value];
	    }
	  else
	    {
	    opacity *= gradient_opacity_constant;
	    }
	}

      // If we have a combined opacity value, then compute the shading
      if ( opacity )
	{
	A_n = *(nptr);
	B_n = *(nptr + Binc);
	C_n = *(nptr + Cinc);
	D_n = *(nptr + Dinc);
	E_n = *(nptr + Einc);
	F_n = *(nptr + Finc);
	G_n = *(nptr + Ginc);
	H_n = *(nptr + Hinc);
	
	final_rd = 
	  red_d_shade[ A_n ] * tA + red_d_shade[ B_n ] * tB + 	
	  red_d_shade[ C_n ] * tC + red_d_shade[ D_n ] * tD + 
	  red_d_shade[ E_n ] * tE + red_d_shade[ F_n ] * tF +	
	  red_d_shade[ G_n ] * tG + red_d_shade[ H_n ] * tH;
	
	final_gd = 
	  green_d_shade[ A_n ] * tA + green_d_shade[ B_n ] * tB + 	
	  green_d_shade[ C_n ] * tC + green_d_shade[ D_n ] * tD + 
	  green_d_shade[ E_n ] * tE + green_d_shade[ F_n ] * tF +	
	  green_d_shade[ G_n ] * tG + green_d_shade[ H_n ] * tH;
	
	final_bd = 
	  blue_d_shade[ A_n ] * tA + blue_d_shade[ B_n ] * tB + 	
	  blue_d_shade[ C_n ] * tC + blue_d_shade[ D_n ] * tD + 
	  blue_d_shade[ E_n ] * tE + blue_d_shade[ F_n ] * tF +	
	  blue_d_shade[ G_n ] * tG + blue_d_shade[ H_n ] * tH;
	
	final_rs = 
	  red_s_shade[ A_n ] * tA + red_s_shade[ B_n ] * tB + 	
	  red_s_shade[ C_n ] * tC + red_s_shade[ D_n ] * tD + 
	  red_s_shade[ E_n ] * tE + red_s_shade[ F_n ] * tF +	
	  red_s_shade[ G_n ] * tG + red_s_shade[ H_n ] * tH;
	
	final_gs = 
	  green_s_shade[ A_n ] * tA + green_s_shade[ B_n ] * tB + 	
	  green_s_shade[ C_n ] * tC + green_s_shade[ D_n ] * tD + 
	  green_s_shade[ E_n ] * tE + green_s_shade[ F_n ] * tF +	
	  green_s_shade[ G_n ] * tG + green_s_shade[ H_n ] * tH;
	
	final_bs = 
	  blue_s_shade[ A_n ] * tA + blue_s_shade[ B_n ] * tB + 	
	  blue_s_shade[ C_n ] * tC + blue_s_shade[ D_n ] * tD + 
	  blue_s_shade[ E_n ] * tE + blue_s_shade[ F_n ] * tF +	
	  blue_s_shade[ G_n ] * tG + blue_s_shade[ H_n ] * tH;
	
	r = CTF[(scalar_value) * 3    ];
	g = CTF[(scalar_value) * 3 + 1];
	b = CTF[(scalar_value) * 3 + 2];
	
	// For this sample we have do not yet have any opacity or
	// shaded intensity yet
	red_shaded_value   = opacity * ( final_rd * r + final_rs );
	green_shaded_value = opacity * ( final_gd * g + final_gs );
	blue_shaded_value  = opacity * ( final_bd * b + final_bs );
	
	// Accumulate intensity and opacity for this sample location   
	accum_red_intensity   += red_shaded_value   * remaining_opacity;
	accum_green_intensity += green_shaded_value * remaining_opacity;
	accum_blue_intensity  += blue_shaded_value  * remaining_opacity;
	remaining_opacity *= (1.0 - opacity);
	}

      // Increment our position and compute our voxel location
      ray_position[0] += ray_increment[0];
      ray_position[1] += ray_increment[1];
      ray_position[2] += ray_increment[2];      
      voxel[0] = (int)( ray_position[0] );
      voxel[1] = (int)( ray_position[1] );
      voxel[2] = (int)( ray_position[2] );
      }
    }

  // Cap the accumulated intensity at 1.0
  if ( accum_red_intensity > 1.0 )
    {
    accum_red_intensity = 1.0;
    }
  if ( accum_green_intensity > 1.0 )
    {
    accum_green_intensity = 1.0;
    }
  if ( accum_blue_intensity > 1.0 )
    {
    accum_blue_intensity = 1.0;
    }
  
  if( remaining_opacity < VTK_REMAINING_OPACITY )
    {
    remaining_opacity = 0.0;
    }

  // Set the return pixel value.  The depth value is currently useless and
  // should be fixed.  What should depth be in this case?  First 
  // non-opaque or total opacity or what??
  rayInfo->RayColor[0] = accum_red_intensity;
  rayInfo->RayColor[1] = accum_green_intensity;
  rayInfo->RayColor[2] = accum_blue_intensity;
  rayInfo->RayColor[3] = 1.0 - remaining_opacity;
  rayInfo->VolumeRayStepsTaken = steps_this_ray;

  if ( remaining_opacity < 1.0 )
    {
    rayInfo->RayDepth = 0.0;
    }
  else 
    {
    rayInfo->RayDepth = VTK_LARGE_FLOAT;
    }
}

// Constructor for the vtkVolumeRayCastCompositeFunction class
vtkVolumeRayCastCompositeFunction::vtkVolumeRayCastCompositeFunction()
{
}

// Destruct the vtkVolumeRayCastCompositeFunction
vtkVolumeRayCastCompositeFunction::~vtkVolumeRayCastCompositeFunction()
{
}

// This is called from RenderAnImage (in vtkDepthPARCMapper.cxx)
// It uses the integer data type flag that is passed in to
// determine what type of ray needs to be cast (which is handled
// by a templated function.  It also uses the shading and
// interpolation types to determine which templated function
// to call.
void vtkVolumeRayCastCompositeFunction::CastRay( struct VolumeRayCastRayInfoStruct *rayInfo,
						 struct VolumeRayCastVolumeInfoStruct *volumeInfo )
{
  void *data_ptr;

  data_ptr = volumeInfo->ScalarDataPointer;

  // Cast the ray for the data type and shading/interpolation type
  if ( volumeInfo->InterpolationType == VTK_NEAREST_INTERPOLATION )
    {
    if ( volumeInfo->Shading == 0 )
      {
      // Nearest neighbor and no shading
      switch ( volumeInfo->ScalarDataType )
	{
	case VTK_UNSIGNED_CHAR:
	  CastRay_NN_Unshaded( this, (unsigned char *)data_ptr, rayInfo, volumeInfo );
	  break;
	case VTK_UNSIGNED_SHORT:
	  CastRay_NN_Unshaded( this, (unsigned short *)data_ptr, rayInfo, volumeInfo );
	  break;
	}
      }
    else
      {
      // Nearest neighbor and shading
      switch ( volumeInfo->ScalarDataType )
	{
	case VTK_UNSIGNED_CHAR:
	  CastRay_NN_Shaded( this, (unsigned char *)data_ptr, rayInfo, volumeInfo );
	  break;
	case VTK_UNSIGNED_SHORT:
	  CastRay_NN_Shaded( this, (unsigned short *)data_ptr, rayInfo, volumeInfo );
	  break;
	}
      }
    }
  else 
    {
    if ( volumeInfo->Shading == 0 )
      {
      // Trilinear interpolation at vertices and no shading
      switch ( volumeInfo->ScalarDataType )
	{
	case VTK_UNSIGNED_CHAR:
	  CastRay_TrilinSample_Unshaded( this, (unsigned char *)data_ptr,  
					 rayInfo, volumeInfo );
	  break;
	case VTK_UNSIGNED_SHORT:
	  CastRay_TrilinSample_Unshaded( this, (unsigned short *)data_ptr, 
					 rayInfo, volumeInfo );
	  break;
	}
      }	
    else
      {
      // Trilinear interpolation and shading
      switch ( volumeInfo->ScalarDataType )
	{
	case VTK_UNSIGNED_CHAR:
	  CastRay_TrilinSample_Shaded( this, (unsigned char *)data_ptr, 
				       rayInfo, volumeInfo );
	  break;
	case VTK_UNSIGNED_SHORT:
	  CastRay_TrilinSample_Shaded( this, (unsigned short *)data_ptr, 
				       rayInfo, volumeInfo );
	  break;
	}
      }	
    }
}

float vtkVolumeRayCastCompositeFunction::GetZeroOpacityThreshold( vtkVolume 
								  *vol )
{
  return vol->GetVolumeProperty()->GetScalarOpacity()->GetFirstNonZeroValue();
}

// We don't need to do any specific initialization here...
void vtkVolumeRayCastCompositeFunction::SpecificFunctionInitialize( 
				vtkRenderer *vtkNotUsed(ren), 
				vtkVolume *vtkNotUsed(vol),
				struct VolumeRayCastVolumeInfoStruct *vtkNotUsed(volumeInfo),
				vtkVolumeRayCastMapper *vtkNotUsed(mapper) )
{
}


// Print method for vtkVolumeRayCastCompositeFunction
// Since there is nothing local to print, just print the object stuff.
void vtkVolumeRayCastCompositeFunction::PrintSelf(ostream& os, vtkIndent indent)
{
  vtkObject::PrintSelf(os,indent);
}





