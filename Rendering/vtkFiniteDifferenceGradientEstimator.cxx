/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkFiniteDifferenceGradientEstimator.cxx
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
#include "vtkFiniteDifferenceGradientEstimator.h"

#include "vtkCharArray.h"
#include "vtkDirectionEncoder.h"
#include "vtkDoubleArray.h"
#include "vtkFloatArray.h"
#include "vtkImageData.h"
#include "vtkIntArray.h"
#include "vtkLongArray.h"
#include "vtkMultiThreader.h"
#include "vtkObjectFactory.h"
#include "vtkShortArray.h"
#include "vtkUnsignedCharArray.h"
#include "vtkUnsignedIntArray.h"
#include "vtkUnsignedLongArray.h"
#include "vtkUnsignedShortArray.h"

#include <math.h>

vtkCxxRevisionMacro(vtkFiniteDifferenceGradientEstimator, "1.34");
vtkStandardNewMacro(vtkFiniteDifferenceGradientEstimator);

// This is the templated function that actually computes the EncodedNormal
// and the GradientMagnitude
template <class T>
void vtkComputeGradients( 
  vtkFiniteDifferenceGradientEstimator *estimator, T *data_ptr,
  int thread_id, int thread_count )
{
  int                 xstep, ystep, zstep;
  int                 x, y, z;
  int                 offset;
  int                 x_start, x_limit;
  int                 y_start, y_limit;
  int                 z_start, z_limit;
  int                 useClip;
  int                 *clip;
  T                   *dptr;
  unsigned char       *gptr;
  unsigned short      *nptr;
  float               n[3], t;
  float               gvalue;
  float               zeroNormalThreshold;
  int                 useBounds;
  int                 bounds[6];
  int                 size[3];
  float               aspect[3];
  int                 xlow, xhigh;
  float               scale, bias;
  int                 computeGradientMagnitudes;
  vtkDirectionEncoder *direction_encoder;
  int                 zeroPad;
  
  estimator->GetInputSize( size );
  estimator->GetInputAspect( aspect );
  computeGradientMagnitudes = estimator->GetComputeGradientMagnitudes();
  scale = estimator->GetGradientMagnitudeScale();
  bias = estimator->GetGradientMagnitudeBias();
  zeroPad = estimator->GetZeroPad();
  
  // adjust the aspect
  aspect[0] = aspect[0] * 2.0 * estimator->SampleSpacingInVoxels;
  aspect[1] = aspect[1] * 2.0 * estimator->SampleSpacingInVoxels;
  aspect[2] = aspect[2] * 2.0 * estimator->SampleSpacingInVoxels;
  
  // Compute steps through the volume in x, y, and z
  xstep = 1;
  ystep = size[0];
  zstep = size[0] * size[1];

  // Multiply by the spacing used for normal estimation
  xstep *= estimator->SampleSpacingInVoxels;
  ystep *= estimator->SampleSpacingInVoxels;
  zstep *= estimator->SampleSpacingInVoxels;

  // Get the length at or below which normals are considered to
  // be "zero"
  zeroNormalThreshold = estimator->GetZeroNormalThreshold();
  
  useBounds = estimator->GetBoundsClip();
  
  // Compute an offset based on the thread_id. The volume will
  // be broken into large slabs (thread_count slabs). For this thread
  // we need to access the correct slab. Also compute the z plane that
  // this slab starts on, and the z limit of this slab (one past the
  // end of the slab)
  if ( useBounds )
    {
    estimator->GetBounds( bounds );
    x_start = bounds[0];
    x_limit = bounds[1]+1;
    y_start = bounds[2];
    y_limit = bounds[3]+1;
    z_start = (int)(( (float)thread_id / (float)thread_count ) *
                    (float)(bounds[5]-bounds[4]+1) ) + bounds[4];
    z_limit = (int)(( (float)(thread_id + 1) / (float)thread_count ) *
                    (float)(bounds[5]-bounds[4]+1) ) + bounds[4];
    }
  else
    {
    x_start = 0;
    x_limit = size[0];
    y_start = 0;
    y_limit = size[1];
    z_start = (int)(( (float)thread_id / (float)thread_count ) *
                    size[2] );
    z_limit = (int)(( (float)(thread_id + 1) / (float)thread_count ) *
                    size[2] );
    }

  // Do final error checking on limits - make sure they are all within bounds
  // of the scalar input
  
  x_start = (x_start<0)?(0):(x_start);
  y_start = (y_start<0)?(0):(y_start);
  z_start = (z_start<0)?(0):(z_start);
  
  x_limit = (x_limit>size[0])?(size[0]):(x_limit);
  y_limit = (y_limit>size[1])?(size[1]):(y_limit);
  z_limit = (z_limit>size[2])?(size[2]):(z_limit);


  direction_encoder = estimator->GetDirectionEncoder();

  useClip = estimator->GetUseCylinderClip();
  clip = estimator->GetCircleLimits();

  // Loop through all the data and compute the encoded normal and
  // gradient magnitude for each scalar location
  for ( z = z_start; z < z_limit; z++ )
    {
    for ( y = y_start; y < y_limit; y++ )
      {
      if ( useClip )
        {
        xlow = ((clip[2*y])>x_start)?(clip[2*y]):(x_start);
        xhigh = ((clip[2*y+1]+1)<x_limit)?(clip[2*y+1]+1):(x_limit);
        }
      else
        {
        xlow = x_start;
        xhigh = x_limit;
        }
      offset = z * size[0] * size[1] + y * size[0] + xlow;
      
      // Set some pointers
      dptr = data_ptr + offset;
      nptr = estimator->EncodedNormals + offset;
      gptr = estimator->GradientMagnitudes + offset;

      for ( x = xlow; x < xhigh; x++ )
        {

        // Use a central difference method if possible,
        // otherwise use a forward or backward difference if
        // we are on the edge

        // Compute the X component
        if ( x < estimator->SampleSpacingInVoxels ) 
          {
          if ( zeroPad )
            {
            n[0] = -((float)*(dptr+xstep));
            }
          else
            {
            n[0] = 2.0*((float)*(dptr) - (float)*(dptr+xstep));
            }
          }
        else if ( x >= size[0] - estimator->SampleSpacingInVoxels )
          {
          if ( zeroPad )
            {
            n[0] =  ((float)*(dptr-xstep));
            }
          else
            {
            n[0] = 2.0*((float)*(dptr-xstep) - (float)*(dptr));
            }
          }
        else
          {
          n[0] = (float)*(dptr-xstep) - (float)*(dptr+xstep); 
          }
        
        // Compute the Y component
        if ( y < estimator->SampleSpacingInVoxels )
          {
          if ( zeroPad )
            {
            n[1] = -((float)*(dptr+ystep));
            }
          else
            {
            n[1] = 2.0*((float)*(dptr) - (float)*(dptr+ystep)); 
            }
          }
        else if ( y >= size[1] - estimator->SampleSpacingInVoxels )
          {
          if ( zeroPad )
            {
            n[1] =  ((float)*(dptr-ystep));
            }
          else
            {
            n[1] = 2.0*((float)*(dptr-ystep) - (float)*(dptr)); 
            }
          }
        else
          {
          n[1] = (float)*(dptr-ystep) - (float)*(dptr+ystep); 
          }
        
        // Compute the Z component
        if ( z < estimator->SampleSpacingInVoxels )
          {
          if ( zeroPad )
            {
            n[2] = -((float)*(dptr+zstep));
            }
          else
            {
            n[2] = 2.0*((float)*(dptr) - (float)*(dptr+zstep)); 
            }
          }
        else if ( z >= size[2] - estimator->SampleSpacingInVoxels )
          {
          if ( zeroPad )
            {
            n[2] =  ((float)*(dptr-zstep));
            }
          else
            {
            n[2] = 2.0*((float)*(dptr-zstep) - (float)*(dptr)); 
            }
          }
        else
          {
          n[2] = (float)*(dptr-zstep) - (float)*(dptr+zstep); 
          }

        // Take care of the aspect ratio of the data
        // Scaling in the vtkVolume is isotropic, so this is the
        // only place we have to worry about non-isotropic scaling.
        n[0] /= aspect[0];
        n[1] /= aspect[1];
        n[2] /= aspect[2];
        
        // Compute the gradient magnitude
        t = sqrt( (double)( n[0]*n[0] + 
                            n[1]*n[1] + 
                            n[2]*n[2] ) );
        
        if ( computeGradientMagnitudes )
          {
          // Encode this into an 8 bit value 
          gvalue = (t + bias) * scale; 
          
          if ( gvalue < 0.0 )
            {
            *gptr = 0;
            }
          else if ( gvalue > 255.0 )
            {
            *gptr = 255;
            }
          else 
            {
            *gptr = (unsigned char) gvalue;
            }
          gptr++;
          }

        // Normalize the gradient direction
        if ( t > zeroNormalThreshold )
          {
          n[0] /= t;
          n[1] /= t;
          n[2] /= t;
          }
        else
          {
          n[0] = n[1] = n[2] = 0.0;
          }

        // Convert the gradient direction into an encoded index value
        *nptr = direction_encoder->GetEncodedDirection( n );
        nptr++;
        dptr++;

        }
      }
    }
}

// Construct a vtkFiniteDifferenceGradientEstimator 
vtkFiniteDifferenceGradientEstimator::vtkFiniteDifferenceGradientEstimator()
{
  this->SampleSpacingInVoxels      = 1;
}

// Destruct a vtkFiniteDifferenceGradientEstimator - free up any memory used
vtkFiniteDifferenceGradientEstimator::~vtkFiniteDifferenceGradientEstimator()
{
}

static VTK_THREAD_RETURN_TYPE vtkSwitchOnDataType( void *arg )
{
  vtkFiniteDifferenceGradientEstimator   *estimator;
  int                                    thread_count;
  int                                    thread_id;
  vtkDataArray                           *scalars;

  thread_id = ((vtkMultiThreader::ThreadInfoStruct *)(arg))->ThreadID;
  thread_count = ((vtkMultiThreader::ThreadInfoStruct *)(arg))->NumberOfThreads;
  estimator = (vtkFiniteDifferenceGradientEstimator *)
    (((vtkMultiThreader::ThreadInfoStruct *)(arg))->UserData);
  scalars = estimator->Input->GetPointData()->GetScalars();

  if (scalars == NULL)
    {
    return VTK_THREAD_RETURN_VALUE;
    }
  
  // Find the data type of the Input and call the correct 
  // templated function to actually compute the normals and magnitudes
  
  switch ( scalars->GetDataType() )
    {
    case VTK_CHAR:
      {
      char *ptr = ((vtkCharArray *) scalars)->GetPointer(0);
      vtkComputeGradients( estimator, ptr, thread_id, thread_count );
      }
    break;
    case VTK_UNSIGNED_CHAR:
      {
      unsigned char *ptr = ((vtkUnsignedCharArray *) 
                            scalars)->GetPointer(0);
      vtkComputeGradients( estimator, ptr, thread_id, thread_count );
      }
    break;
    case VTK_SHORT:
      {
      short *ptr = ((vtkShortArray *) scalars)->GetPointer(0);
      vtkComputeGradients( estimator, ptr, thread_id, thread_count );
      }
    break;
    case VTK_UNSIGNED_SHORT:
      {
      unsigned short *ptr = ((vtkUnsignedShortArray *) 
                             scalars)->GetPointer(0);
      vtkComputeGradients( estimator, ptr, thread_id, thread_count );
      }
    break;
    case VTK_INT:
      {
      int *ptr = ((vtkIntArray *) scalars)->GetPointer(0);
      vtkComputeGradients( estimator, ptr, thread_id, thread_count );
      }
    break;
    case VTK_UNSIGNED_INT:
      {
      unsigned int *ptr = ((vtkUnsignedIntArray *) 
                           scalars)->GetPointer(0);
      vtkComputeGradients( estimator, ptr, thread_id, thread_count );
      }
    break;
    case VTK_LONG:
      {
      long *ptr = ((vtkLongArray *) scalars)->GetPointer(0);
      vtkComputeGradients( estimator, ptr, thread_id, thread_count );
      }
    break;
    case VTK_UNSIGNED_LONG:
      {
      unsigned long *ptr = ((vtkUnsignedLongArray *) 
                            scalars)->GetPointer(0);
      vtkComputeGradients( estimator, ptr, thread_id, thread_count );
      }
    break;
    case VTK_FLOAT:
      {
      float *ptr = ((vtkFloatArray *) scalars)->GetPointer(0);
      vtkComputeGradients( estimator, ptr, thread_id, thread_count );
      }
    break;
    case VTK_DOUBLE:
      {
      double *ptr = ((vtkDoubleArray *) scalars)->GetPointer(0);
      vtkComputeGradients( estimator, ptr, thread_id, thread_count );
      }
    break;
    default:
      vtkGenericWarningMacro("unable to encode scalar type!");
    }
  
  return VTK_THREAD_RETURN_VALUE;
}


// This method is used to compute the encoded normal and the
// magnitude of the gradient for each voxel location in the 
// Input.
void vtkFiniteDifferenceGradientEstimator::UpdateNormals( )
{
  vtkDebugMacro( << "Updating Normals!" );
  this->Threader->SetNumberOfThreads( this->NumberOfThreads );
  
  this->Threader->SetSingleMethod( vtkSwitchOnDataType,
                                  (vtkObject *)this );
  
  this->Threader->SingleMethodExecute();
}

// Print the vtkFiniteDifferenceGradientEstimator
void vtkFiniteDifferenceGradientEstimator::PrintSelf(ostream& os, 
                                                     vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
  
  os << indent << "Sample spacing in voxels: " << 
    this->SampleSpacingInVoxels << endl;
}
