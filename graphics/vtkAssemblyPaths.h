/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkAssemblyPaths.h
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
// .NAME vtkAssemblyPaths - a list of lists of actors representing an assembly hierarchy
// .SECTION Description
// vtkAssemblyPaths represents a hierarchy of assemblies as a sequence of
// paths. Each path is a list of actors, starting from the root of the
// assembly down to the leaf actors. Methods are also provided to manipulate
// the path including propagating transformation matrices and actor properties.

// .SECTION see also
// vtkAssembly vtkActor

#ifndef __vtkAssemblyPaths_h
#define __vtkAssemblyPaths_h

#include "vtkActorCollection.h"
class vtkActor;

class VTK_EXPORT vtkAssemblyPaths : public vtkCollection
{
public:
  static vtkAssemblyPaths *New();
  vtkTypeMacro(vtkAssemblyPaths,vtkCollection);

  // Description:
  // Add a path to the list.
  void AddItem(vtkActorCollection *a);

  // Description:
  // Remove a path from the list.
  void RemoveItem(vtkActorCollection *a);

  // Description:
  // Determine whether a particular path is present. Returns its position
  // in the list.
  int IsItemPresent(vtkActorCollection *a);

  // Description:
  // Get the next path in the list.
  vtkActorCollection *GetNextItem();

protected:
  vtkAssemblyPaths() {};
  ~vtkAssemblyPaths() {};
  vtkAssemblyPaths(const vtkAssemblyPaths&) {};
  void operator=(const vtkAssemblyPaths&) {};
  
private:
  // hide the standard AddItem from the user and the compiler.
  void AddItem(vtkObject *o) { this->vtkCollection::AddItem(o); };
  void RemoveItem(vtkObject *o) { this->vtkCollection::RemoveItem(o); };
  void RemoveItem(int i) { this->vtkCollection::RemoveItem(i); };
  int  IsItemPresent(vtkObject *o) { return this->vtkCollection::IsItemPresent(o);};
};

inline void vtkAssemblyPaths::AddItem(vtkActorCollection *a) 
{
  this->vtkCollection::AddItem((vtkObject *)a);
}

inline void vtkAssemblyPaths::RemoveItem(vtkActorCollection *a) 
{
  this->vtkCollection::RemoveItem((vtkObject *)a);
}

inline int vtkAssemblyPaths::IsItemPresent(vtkActorCollection *a) 
{
  return this->vtkCollection::IsItemPresent((vtkObject *)a);
}

inline vtkActorCollection *vtkAssemblyPaths::GetNextItem() 
{ 
  return (vtkActorCollection *)(this->GetNextItemAsObject());
}

#endif
