/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkObjectBase.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkObjectBase.h"
#include "vtkDebugLeaks.h"

#define vtkBaseDebugMacro(x)

// avoid dll boundary problems
#ifdef _WIN32
void* vtkObjectBase::operator new(size_t nSize)
{
  void* p=malloc(nSize);
  return p;
}

void vtkObjectBase::operator delete( void *p )
{
  free(p);
}
#endif 

// ------------------------------------vtkObjectBase----------------------
// This operator allows all subclasses of vtkObjectBase to be printed via <<.
// It in turn invokes the Print method, which in turn will invoke the
// PrintSelf method that all objects should define, if they have anything
// interesting to print out.
ostream& operator<<(ostream& os, vtkObjectBase& o)
{
  o.Print(os);
  return os;
}

// Create an object with Debug turned off and modified time initialized 
// to zero.
vtkObjectBase::vtkObjectBase()
{
  this->ReferenceCount = 1;
  // initial reference count = 1 and reference counting on.
}

vtkObjectBase::~vtkObjectBase() 
{
  // warn user if reference counting is on and the object is being referenced
  // by another object
  if ( this->ReferenceCount > 0)
    {
    vtkGenericWarningMacro(<< "Trying to delete object with non-zero reference count.");
    }
}

int vtkObjectBase::IsTypeOf(const char *name) 
{
  if ( !strcmp("vtkObjectBase",name) )
    {
    return 1;
    }
  return 0;
}

int vtkObjectBase::IsA(const char *type)
{
  return this->vtkObjectBase::IsTypeOf(type);
}

// Delete a vtk object. This method should always be used to delete an object 
// when the new operator was used to create it. Using the C++ delete method
// will not work with reference counting.
void vtkObjectBase::Delete() 
{
  this->UnRegister((vtkObjectBase *)NULL);
}

void vtkObjectBase::Print(ostream& os)
{
  vtkIndent indent;

  this->PrintHeader(os,0); 
  this->PrintSelf(os, indent.GetNextIndent());
  this->PrintTrailer(os,0);
}

void vtkObjectBase::PrintHeader(ostream& os, vtkIndent indent)
{
  os << indent << this->GetClassName() << " (" << this << ")\n";
}

// Chaining method to print an object's instance variables, as well as
// its superclasses.
void vtkObjectBase::PrintSelf(ostream& os, vtkIndent indent)
{
  os << indent << "Reference Count: " << this->ReferenceCount << "\n";
}

void vtkObjectBase::PrintTrailer(ostream& os, vtkIndent indent)
{
  os << indent << "\n";
}

// Description:
// Sets the reference count (use with care)
void vtkObjectBase::SetReferenceCount(int ref)
{
  this->ReferenceCount = ref;
  vtkBaseDebugMacro(<< "Reference Count set to " << this->ReferenceCount);
}

// Description:
// Increase the reference count (mark as used by another object).
void vtkObjectBase::Register(vtkObjectBase*)
{
  this->ReferenceCount++;
  if (this->ReferenceCount <= 0)
    {
    delete this;
    }
}

// Description:
// Decrease the reference count (release by another object).
void vtkObjectBase::UnRegister(vtkObjectBase* o)
{
  if (o)
    {
    vtkBaseDebugMacro(
      << "UnRegistered by "
      << o->GetClassName() << " (" << o << "), ReferenceCount = "
      << (this->ReferenceCount-1));
    }
  else
    {
    vtkBaseDebugMacro(
      << "UnRegistered " << this->GetClassName() 
      << " by NULL, ReferenceCount = "
      << (this->ReferenceCount-1));
    }

  if (--this->ReferenceCount <= 0)
    {
#ifdef VTK_DEBUG_LEAKS
    vtkDebugLeaks::DestructClass(this->GetClassName());
#endif
    // invoke the delete method
    // Here we should call delete method
    delete this;
    }
}

void vtkObjectBase::CollectRevisions(ostream& os)
{
  os << "vtkObjectBase 1.7\n";
}

void vtkObjectBase::PrintRevisions(ostream& os)
{
  ostrstream revisions;
  this->CollectRevisions(revisions);
  revisions << ends;
  const char* c = revisions.str();
  while(*c)
    {
    const char* beginClass = 0;
    const char* endClass = 0;
    const char* beginRevision = 0;
    const char* endRevision = 0;
    for(;*c && *c != '\n'; ++c)
      {
      if(!beginClass && *c != ' ')
        {
        beginClass = c;
        }
      else if(beginClass && !endClass && *c == ' ')
        {
        endClass = c;
        }
      else if(endClass && !beginRevision && (*c >= '0' && *c <= '9'))
        {
        beginRevision = c;
        }
      else if(beginRevision && !endRevision && *c == ' ')
        {
        endRevision = c;
        }
      }
    if(beginClass && endClass && beginRevision && endRevision)
      {
      os.write(beginClass, endClass-beginClass);
      os << " ";
      os.write(beginRevision, endRevision-beginRevision);
      os << "\n";
      }
    if(*c == '\n')
      {
      ++c;
      }
    }
  revisions.rdbuf()->freeze(0);
}

//----------------------------------------------------------------------------
void vtkObjectBase::ReportReferences(vtkGarbageCollector*)
{
  // vtkObjectBase has no references to report.
}

//----------------------------------------------------------------------------
void vtkObjectBase::RemoveReferences()
{
  // vtkObjectBase has no references to remove.
}

//----------------------------------------------------------------------------
void vtkObjectBase::GarbageCollectionStarting()
{
  // Do not delete this object until garbage collection is finishing.
  this->Register(0);
}

//----------------------------------------------------------------------------
void vtkObjectBase::GarbageCollectionFinishing()
{
  // Delete this object now.
  this->UnRegister(0);
}
