/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkFloatVectors.h
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
// .NAME vtkFloatVectors - (obsolete)floating point representation of 3D vectors
// .SECTION Description
// vtkFloatVectors is a concrete implementation of vtkVectors. Vectors are
// represented using float values.

#ifndef __vtkFloatVectors_h
#define __vtkFloatVectors_h

#include "vtkVectors.h"
#include "vtkFloatArray.h"

class VTK_EXPORT vtkFloatVectors : public vtkVectors
{
public:
  static vtkFloatVectors *New();
  
  // Description:
  // Set the data type for this object.
  void SetDataType(int dataType);

  // Description:
  // Set the data for this object. Only accepts VTK_FLOAT type.
  void SetData(vtkDataArray *);

  // Description:
  // Get pointer to array of data starting at data position "id".
  float *GetPointer(const int id);

  // Description:
  // Get pointer to data array. Useful for direct writes of data. MaxId is
  // bumped by number (and memory allocated if necessary). Id is the
  // location you wish to write into; number is the number of vectors to
  // write.
  float *WritePointer(const int id, const int number);

protected:
  vtkFloatVectors() {};
  ~vtkFloatVectors() {};
  vtkFloatVectors(const vtkFloatVectors&) {};
  void operator=(const vtkFloatVectors&) {};
  
};

inline float *vtkFloatVectors::GetPointer(const int id)
{
  return ((vtkFloatArray *)this->Data)->GetPointer(3*id);
} 

inline float *vtkFloatVectors::WritePointer(const int id, const int number)
{
  return ((vtkFloatArray *)this->Data)->WritePointer(3*id,3*number);
}

inline void vtkFloatVectors::SetData(vtkDataArray *data)
{
  if ( data->GetDataType() != VTK_FLOAT )
    {
    vtkErrorMacro(<<"Float vectors only accepts float data type");
    return;
    }
  vtkVectors::SetData(data);
}

inline void vtkFloatVectors::SetDataType(int type)
{
  if ( type != VTK_FLOAT )
    {
    vtkErrorMacro(<<"Float vectors only accepts float data type");
    return;
    }

  vtkVectors::SetDataType(type);
}


#endif
