/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkCharArray.h
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
// .NAME vtkCharArray - dynamic, self-adjusting character array
// .SECTION Description
// vtkCharArray is an array of character values. It provides 
// methods for insertion and retrieval of characters, and will automatically 
// resize itself to hold new data.

#ifndef __vtkCharArray_h
#define __vtkCharArray_h

#include "vtkDataArray.h"

class VTK_EXPORT vtkCharArray : public vtkDataArray
{
public:
  vtkCharArray(int numComp=1);
  ~vtkCharArray();
  int Allocate(const int sz, const int ext=1000);
  void Initialize();
  static vtkCharArray *New() {return new vtkCharArray;};
  const char *GetClassName() {return "vtkCharArray";};
  void PrintSelf(ostream& os, vtkIndent indent);

  // satisfy vtkDataArray API
  vtkDataArray *MakeObject() {return new vtkCharArray(this->NumberOfComponents);};
  int GetDataType() {return VTK_CHAR;};
  void SetNumberOfTuples(const int number);
  float *GetTuple(const int i);
  void GetTuple(const int i, float tuple[]);
  void SetTuple(const int i, const float tuple[]);
  void InsertTuple(const int i, const float tuple[]);
  int InsertNextTuple(const float tuple[]);
  void Squeeze();

  // overload vtkDataArray for efficiency
  float GetComponent(const int i, const int j);
  void SetComponent(const int i, const int j, const float c);
  void InsertComponent(const int i, const int j, const float c);

  // native access/insertion methods
  char GetValue(const int id);
  void SetNumberOfValues(const int number);
  void SetValue(const int id, const char c);
  void InsertValue(const int id, const char c);
  int InsertNextValue(const char c);
  char *GetPointer(const int id) {return this->Array + id;}

  char *WritePointer(const int id, const int number);
  void *GetVoidPointer(const int id) {return (void *)this->GetPointer(id);};
  void DeepCopy(vtkCharArray& ia);

private:
  char *Array;    // pointer to data
  char *Resize(const int sz);  // function to resize data

  int TupleSize; //used for data conversion
  float *Tuple;
};

// Description:
// Get the data at a particular index.
inline char vtkCharArray::GetValue(const int id) {return this->Array[id];}

// Description:
// Specify the number of values for this object to hold. Does an
// allocation as well as setting the MaxId ivar. Used in conjunction with
// SetValue() method for fast insertion.
inline void vtkCharArray::SetNumberOfValues(const int number) 
{
  this->Allocate(number);
  this->MaxId = number - 1;
}

// Description:
// Set the data at a particular index. Does not do range checking. Make sure
// you use the method SetNumberOfValues() before inserting data.
inline void vtkCharArray::SetValue(const int id, const char value) 
{
  this->Array[id] = value;
}

// Description:
// Get the address of a particular data index. Make sure data is allocated
// for the number of items requested. Set MaxId according to the number of
// data values requested.
inline char *vtkCharArray::WritePointer(const int id, const int number) 
{
  int newSize=id+number;
  if ( newSize > this->Size ) this->Resize(newSize);
  if ( (--newSize) > this->MaxId ) this->MaxId = newSize;
  return this->Array + id;
}

// Description:
// Insert data at a specified position in the array.
inline void vtkCharArray::InsertValue(const int id, const char c)
{
  if ( id >= this->Size ) this->Resize(id+1);
  this->Array[id] = c;
  if ( id > this->MaxId ) this->MaxId = id;
}

// Description:
// Insert data at the end of the array. Return its location in the array.
inline int vtkCharArray::InsertNextValue(const char c)
{
  this->InsertValue (++this->MaxId,c); 
  return this->MaxId;
}

// Description:
// Resize object to just fit data requirement. Reclaims extra memory.
inline void vtkCharArray::Squeeze() {this->Resize (this->MaxId+1);}

#endif
