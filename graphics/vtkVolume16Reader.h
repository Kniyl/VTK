/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkVolume16Reader.h
  Language:  C++
  Date:      $Date$
  Version:   $Revision$


Copyright (c) 1993-2001 Ken Martin, Will Schroeder, Bill Lorensen 
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.

 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

 * Neither name of Ken Martin, Will Schroeder, or Bill Lorensen nor the names
   of any contributors may be used to endorse or promote products derived
   from this software without specific prior written permission.

 * Modified source versions must be plainly marked as such, and must not be
   misrepresented as being the original software.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS ``AS IS''
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=========================================================================*/
// .NAME vtkVolume16Reader - read 16 bit image files
// .SECTION Description
// vtkVolume16Reader is a source object that reads 16 bit image files.
//
// Volume16Reader creates structured point datasets. The dimension of the 
// dataset depends upon the number of files read. Reading a single file 
// results in a 2D image, while reading more than one file results in a 
// 3D volume.
//
// File names are created using FilePattern and FilePrefix as follows:
// sprintf (filename, FilePattern, FilePrefix, number);
// where number is in the range ImageRange[0] to ImageRange[1]. If
// ImageRange[1] <= ImageRange[0], then slice number ImageRange[0] is
// read. Thus to read an image set ImageRange[0] = ImageRange[1] = slice 
// number. The default behavior is to read a single file (i.e., image slice 1).
//
// The DataMask instance variable is used to read data files with imbedded
// connectivity or segmentation information. For example, some data has
// the high order bit set to indicate connected surface. The DataMask allows
// you to select this data. Other important ivars include HeaderSize, which
// allows you to skip over initial info, and SwapBytes, which turns on/off
// byte swapping.
//
// The Transform instance variable specifies a permutation transformation
// to map slice space into world space. vtkImageReader has replaced the
// functionality of this class and should be used instead.

// .SECTION See Also
// vtkSliceCubes vtkMarchingCubes vtkImageReader

#ifndef __vtkVolume16Reader_h
#define __vtkVolume16Reader_h

#include <stdio.h>
#include "vtkVolumeReader.h"
#include "vtkTransform.h"

#define VTK_FILE_BYTE_ORDER_BIG_ENDIAN 0
#define VTK_FILE_BYTE_ORDER_LITTLE_ENDIAN 1

class VTK_EXPORT vtkVolume16Reader : public vtkVolumeReader
{
public:
  vtkTypeMacro(vtkVolume16Reader,vtkVolumeReader);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Construct object with NULL file prefix; file pattern "%s.%d"; image range 
  // set to (1,1); data origin (0,0,0); data spacing (1,1,1); no data mask;
  // header size 0; and byte swapping turned off.
  static vtkVolume16Reader *New();

  // Description:
  // Specify the dimensions for the data.
  vtkSetVector2Macro(DataDimensions,int);
  vtkGetVectorMacro(DataDimensions,int,2);

  // Description:
  // Specify a mask used to eliminate data in the data file (e.g.,
  // connectivity bits).
  vtkSetMacro(DataMask,unsigned short);
  vtkGetMacro(DataMask,unsigned short);

  // Description:
  // Specify the number of bytes to seek over at start of image.
  vtkSetMacro(HeaderSize,int);
  vtkGetMacro(HeaderSize,int);

  // Description:
  // These methods should be used instead of the SwapBytes methods.
  // They indicate the byte ordering of the file you are trying
  // to read in. These methods will then either swap or not swap
  // the bytes depending on the byte ordering of the machine it is
  // being run on. For example, reading in a BigEndian file on a
  // BigEndian machine will result in no swapping. Trying to read
  // the same file on a LittleEndian machine will result in swapping.
  // As a quick note most UNIX machines are BigEndian while PC's
  // and VAX tend to be LittleEndian. So if the file you are reading
  // in was generated on a VAX or PC, SetDataByteOrderToLittleEndian otherwise
  // SetDataByteOrderToBigEndian. 
  void SetDataByteOrderToBigEndian();
  void SetDataByteOrderToLittleEndian();
  int GetDataByteOrder();
  void SetDataByteOrder(int);
  const char *GetDataByteOrderAsString();

  // Description:
  // Turn on/off byte swapping.
  vtkSetMacro(SwapBytes,int);
  vtkGetMacro(SwapBytes,int);
  vtkBooleanMacro(SwapBytes,int);

  // Description:
  // Set/Get transformation matrix to transform the data from slice space
  // into world space. This matrix must be a permutation matrix. To qualify,
  // the sums of the rows must be + or - 1.
  vtkSetObjectMacro(Transform,vtkTransform);
  vtkGetObjectMacro(Transform,vtkTransform);

  // Description:
  // Other objects make use of these methods
  vtkStructuredPoints *GetImage(int ImageNumber);

protected:
  vtkVolume16Reader();
  ~vtkVolume16Reader();
  vtkVolume16Reader(const vtkVolume16Reader&) {};
  void operator=(const vtkVolume16Reader&) {};

  void Execute();
  void ExecuteInformation();
  int   DataDimensions[2];
  unsigned short DataMask;
  int   SwapBytes;
  int   HeaderSize;
  vtkTransform *Transform;

  void TransformSlice (unsigned short *slice, unsigned short *pixels, int k, int dimensions[3], int bounds[3]);
  void ComputeTransformedDimensions(int dimensions[3]);
  void ComputeTransformedBounds(int bounds[6]);
  void ComputeTransformedSpacing(float Spacing[3]);
  void ComputeTransformedOrigin(float origin[3]);
  void AdjustSpacingAndOrigin(int dimensions[3], float Spacing[3], float origin[3]);
  vtkScalars *ReadImage(int ImageNumber);
  vtkScalars *ReadVolume(int FirstImage, int LastImage);
  int Read16BitImage(FILE *fp, unsigned short *pixels, int xsize, int ysize, 
		     int skip, int swapBytes);

};

#endif


