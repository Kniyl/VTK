/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkImagingFactory.cxx
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

#include "vtkObjectFactory.h"
#include "vtkImagingFactory.h"
#include "vtkToolkits.h"
#include "vtkDebugLeaks.h"

#ifdef VTK_USE_OGLR
#include "vtkOpenGLImageMapper.h"
#include "vtkXOpenGLTextMapper.h"
#include "vtkOpenGLPolyDataMapper2D.h"
#endif

#if defined(VTK_MANGLE_MESA)
#include "vtkMesaImageMapper.h"
#include "vtkXMesaTextMapper.h"
#include "vtkMesaPolyDataMapper2D.h"
#endif

#ifdef _WIN32
#include "vtkOpenGLImageMapper.h"
#include "vtkWin32OpenGLTextMapper.h"
#include "vtkOpenGLPolyDataMapper2D.h"
#endif

#ifdef VTK_USE_CARBON
#include "vtkOpenGLImageMapper.h"
#include "vtkCarbonTextMapper.h"
#include "vtkOpenGLPolyDataMapper2D.h"
#endif
#ifdef VTK_USE_COCOA
#include "vtkOpenGLImageMapper.h"
#include "vtkCocoaTextMapper.h"
#include "vtkOpenGLPolyDataMapper2D.h"
#endif

#include "vtkCriticalSection.h"
static vtkSimpleCriticalSection vtkUseMesaClassesCriticalSection;
int vtkImagingFactory::UseMesaClasses = 0;

vtkCxxRevisionMacro(vtkImagingFactory, "1.23");


const char *vtkImagingFactoryGetRenderLibrary()
{
  const char *temp;
  
  // first check the environment variable
  temp = getenv("VTK_RENDERER");
  
  // Backward compatibility
  if ( temp )
    {
    if (!strcmp("oglr",temp))
      {
      temp = "OpenGL";
      }
    else if (!strcmp("woglr",temp))
      {
      temp = "Win32OpenGL";
      }
    else if (strcmp("Mesa",temp) && strcmp("OpenGL",temp) && 
             strcmp("Win32OpenGL",temp))
      {
      vtkGenericWarningMacro(<<"VTK_RENDERER set to unsupported type:" << temp);
      temp = NULL;
      }
    }
  
  // if nothing is set then work down the list of possible renderers
  if ( !temp )
    {
#ifdef VTK_USE_OGLR
    temp = "OpenGL";
#endif
#ifdef _WIN32
    temp = "Win32OpenGL";
#endif
#ifdef VTK_USE_CARBON
    temp = "CarbonOpenGL";
#endif
#ifdef VTK_USE_COCOA
    temp = "CocoaOpenGL";
#endif
    }
  
  return temp;
}

vtkObject* vtkImagingFactory::CreateInstance(const char* vtkclassname )
{
  // first check the object factory
  vtkObject *ret = vtkObjectFactory::CreateInstance(vtkclassname);
  if (ret)
    {
    return ret;
    }
  // if the factory failed to create the object,
  // then destroy it now, as vtkDebugLeaks::ConstructClass was called
  // with vtkclassname, and not the real name of the class
#ifdef VTK_DEBUG_LEAKS
  vtkDebugLeaks::DestructClass(vtkclassname);
#endif

  const char *rl = vtkImagingFactoryGetRenderLibrary();

#ifdef VTK_USE_OGLR
  if (!strcmp("OpenGL",rl))
    {
    if(strcmp(vtkclassname, "vtkTextMapper") == 0)
      {
#if defined(VTK_MANGLE_MESA)
      if ( vtkImagingFactory::UseMesaClasses )
	{
	return vtkXMesaTextMapper::New();
	}
#endif
      return vtkXOpenGLTextMapper::New();
      }
    if(strcmp(vtkclassname, "vtkImageMapper") == 0)
      {
#if defined(VTK_MANGLE_MESA)
      if ( vtkImagingFactory::UseMesaClasses )
	{
	return vtkMesaImageMapper::New();
	}
#endif
      return vtkOpenGLImageMapper::New();
      }
    if(strcmp(vtkclassname, "vtkPolyDataMapper2D") == 0)
      {
#if defined(VTK_MANGLE_MESA)
      if ( vtkImagingFactory::UseMesaClasses )
	{
	return vtkMesaPolyDataMapper2D::New();
	}
#endif
      return vtkOpenGLPolyDataMapper2D::New();
      }
    }
#endif

#ifdef _WIN32
  if (!strcmp("Win32OpenGL",rl))
    {
    if(strcmp(vtkclassname, "vtkTextMapper") == 0)
      {
      return vtkWin32OpenGLTextMapper::New();
      }
    if(strcmp(vtkclassname, "vtkImageMapper") == 0)
      {
      return vtkOpenGLImageMapper::New();
      }
    if(strcmp(vtkclassname, "vtkPolyDataMapper2D") == 0)
      {
      return vtkOpenGLPolyDataMapper2D::New();
      }
    }
#endif

#ifdef VTK_USE_CARBON
  if (!strcmp("CarbonOpenGL",rl))
    {
    if(strcmp(vtkclassname, "vtkTextMapper") == 0)
      {
      return vtkCarbonTextMapper::New();
      }
    if(strcmp(vtkclassname, "vtkImageMapper") == 0)
      {
      return vtkOpenGLImageMapper::New();
      }
    if(strcmp(vtkclassname, "vtkPolyDataMapper2D") == 0)
      {
      return vtkOpenGLPolyDataMapper2D::New();
      }
    }
#endif
#ifdef VTK_USE_COCOA
  if (!strcmp("CocoaOpenGL",rl))
    {
    if(strcmp(vtkclassname, "vtkTextMapper") == 0)
      {
      return vtkCocoaTextMapper::New();
      }
    if(strcmp(vtkclassname, "vtkImageMapper") == 0)
      {
      return vtkOpenGLImageMapper::New();
      }
    if(strcmp(vtkclassname, "vtkPolyDataMapper2D") == 0)
      {
      return vtkOpenGLPolyDataMapper2D::New();
      }
    }
#endif

  return 0;
}

void vtkImagingFactory::SetUseMesaClasses(int use)
{
  vtkUseMesaClassesCriticalSection.Lock();
  vtkImagingFactory::UseMesaClasses = use;
  vtkUseMesaClassesCriticalSection.Unlock();
}

int vtkImagingFactory::GetUseMesaClasses()
{
  return vtkImagingFactory::UseMesaClasses;
}
