/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkArrayMapIterator.h
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
// .NAME vtkArrayMapIterator - an array map iterator

#ifndef __vtkArrayMapIterator_h
#define __vtkArrayMapIterator_h

#include "vtkAbstractIterator.h"

template <class KeyType,class DataType>
class vtkArrayMapIterator : public vtkAbstractIterator<KeyType,DataType>
{
  friend class vtkArrayMap<KeyType,DataType>;

public:
  virtual const char* GetClassName() const { return "vtkArrayMapIterator"; }

  // Description:
  // Retrieve the index of the element.
  // This method returns VTK_OK if key was retrieved correctly.
  int GetKey(KeyType&);

  // Description:
  // Retrieve the data from the iterator. 
  // This method returns VTK_OK if key was retrieved correctly.
  int GetData(DataType&);

  // Description:
  // Initialize the traversal of the container. 
  // Set the iterator to the "beginning" of the container.
  void InitTraversal();

  // Description:
  // Check if the iterator is at the end of the container. Returns 1 for yes
  // and 0 for no.
  int IsDoneWithTraversal();

  // Description:
  // Increment the iterator to the next location.
  // Return VTK_OK if everything is ok.
  int GoToNextItem();

  // Description:
  // Decrement the iterator to the next location.
  // Return VTK_OK if everything is ok.
  int GoToPreviousItem();

  // Description:
  // Go to the "first" item of the map.
  // Return VTK_OK if everything is ok.
  int GoToFirstItem();

  // Description:
  // Go to the "last" item of the map.
  // Return VTK_OK if everything is ok.
  int GoToLastItem();

protected:
  static vtkArrayMapIterator<KeyType,DataType> *New(); 

  vtkArrayMapIterator() {
    this->Index = 0; 
  }
  virtual ~vtkArrayMapIterator() {}

  vtkIdType Index;

private:
  vtkArrayMapIterator(const vtkArrayMapIterator<KeyType,DataType>&); // Not implemented
  void operator=(const vtkArrayMapIterator<KeyType,DataType>&); // Not implemented
};

#ifdef VTK_NO_EXPLICIT_TEMPLATE_INSTANTIATION
#include "vtkArrayMapIterator.txx"
#endif 

#endif
