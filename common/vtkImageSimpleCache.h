/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkImageSimpleCache.h
  Language:  C++
  Date:      $Date$
  Version:   $Revision$
  Thanks:    Thanks to C. Charles Law who developed this class.

Copyright (c) 1993-1995 Ken Martin, Will Schroeder, Bill Lorensen.

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
// .NAME vtkImageSimpleCache - Caches the last region generated
// .SECTION Description
// vtkImageSimpleCache saves the last generated region.
// If a subsequent region is contained in the cached data, the
// cached data is returned with no call to the filters Update method.
// If the new region is not completely contained in the cached data,
// the cached data is not used.


#ifndef __vtkImageSimpleCache_h
#define __vtkImageSimpleCache_h

#include "vtkImageCache.h"

class VTK_EXPORT vtkImageSimpleCache : public vtkImageCache
{
public:
  static vtkImageSimpleCache *New() {return new vtkImageSimpleCache;};

  vtkTypeMacro(vtkImageSimpleCache,vtkImageCache);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // This method updates the region specified by "UpdateExtent".  
  void Update();

  // Description:
  // This is the most common way to obtain data from a cache.
  // After setting the update extent invoke this method and it
  // will return an ImageData instance containing the requested data.
  vtkImageData *UpdateAndReturnData();

  // Description:
  // This method deletes any data in the cache.
  void ReleaseData();

  // Description:
  // Allocates the scalar data required for the current update extent.
  void AllocateData();

  // Description:
  // return the un filled data of the UpdateExtent in this cache.
  vtkImageData *GetData(); 

  // Description:
  // Convenience method to get the range of the scalar data in the
  // current "UpdateExtent". Returns the (min/max) range.  The components
  // are lumped into one range.  If there are no scalars the method will 
  // return (0,1). Note: Update needs to be called first to create the scalars.
  void GetScalarRange(float range[2]);  
  
protected:
  vtkImageSimpleCache();
  ~vtkImageSimpleCache();
  vtkImageSimpleCache(const vtkImageSimpleCache&) {};
  void operator=(const vtkImageSimpleCache&) {};

  vtkImageData *CachedData;
  vtkTimeStamp GenerateTime;
};

#endif


