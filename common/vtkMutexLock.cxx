/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkMutexLock.cxx
  Language:  C++
  Date:      $Date$
  Version:   $Revision$


Copyright (c) 1993-2000 Ken Martin, Will Schroeder, Bill Lorensen.

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
#include "vtkMutexLock.h"
#include "vtkObjectFactory.h"



//------------------------------------------------------------------------------
vtkMutexLock* vtkMutexLock::New()
{
  // First try to create the object from the vtkObjectFactory
  vtkObject* ret = vtkObjectFactory::CreateInstance("vtkMutexLock");
  if(ret)
    {
    return (vtkMutexLock*)ret;
    }
  // If the factory was unable to create the object, then create it here.
  return new vtkMutexLock;
}




// Construct a new vtkMutexLock 
vtkSimpleMutexLock::vtkSimpleMutexLock()
{
#ifdef VTK_USE_SPROC
  init_lock( &this->MutexLock );
#endif

#ifdef _WIN32
  this->MutexLock = CreateMutex( NULL, FALSE, NULL ); 
#endif

#ifdef VTK_USE_PTHREADS
#ifdef VTK_HP_PTHREADS
  pthread_mutex_init(&(this->MutexLock), pthread_mutexattr_default);
#else
  pthread_mutex_init(&(this->MutexLock), NULL);
#endif
#endif

}

// Destruct the vtkMutexVariable
vtkSimpleMutexLock::~vtkSimpleMutexLock()
{
#ifdef _WIN32
  CloseHandle(this->MutexLock);
#endif

#ifdef VTK_USE_PTHREADS
  pthread_mutex_destroy( &this->MutexLock);
#endif
}

// Lock the vtkMutexLock
void vtkSimpleMutexLock::Lock()
{
#ifdef VTK_USE_SPROC
  spin_lock( &this->MutexLock );
#endif

#ifdef _WIN32
  WaitForSingleObject( this->MutexLock, INFINITE );
#endif

#ifdef VTK_USE_PTHREADS
  pthread_mutex_lock( &this->MutexLock);
#endif
}

// Unlock the vtkMutexLock
void vtkSimpleMutexLock::Unlock()
{
#ifdef VTK_USE_SPROC
  release_lock( &this->MutexLock );
#endif

#ifdef _WIN32
  ReleaseMutex( this->MutexLock );
#endif

#ifdef VTK_USE_PTHREADS
  pthread_mutex_unlock( &this->MutexLock);
#endif
}

int vtkSimpleMutexLock::IsA(const char *type)
{
  if ( !strcmp(this->vtkSimpleMutexLock::GetClassName(),type) )
    {
    return 1;
    }
  return 0;
}

vtkSimpleMutexLock *vtkSimpleMutexLock::SafeDownCast(vtkSimpleMutexLock *o)
{
  return (vtkSimpleMutexLock *)o;
}
  
void vtkMutexLock::PrintSelf(ostream& os, vtkIndent indent)
{
  vtkObject::PrintSelf(os, indent);
}

