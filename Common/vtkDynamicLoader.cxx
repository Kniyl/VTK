/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkDynamicLoader.cxx
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
#include "vtkDynamicLoader.h"

// This file is actually 3 different implementations.
// 1. HP machines which uses shl_load
// 2. Apple OSX which uses NSLinkModule
// 3. Windows which uses LoadLibrary
// 4. Most unix systems which use dlopen (default )
// Each part of the ifdef contains a complete implementation for
// the static methods of vtkDynamicLoader.  

vtkCxxRevisionMacro(vtkDynamicLoader, "1.9");

// ---------------------------------------------------------------
// 1. Implementation for HPUX  machines
#ifdef __hpux
#define VTKDYNAMICLOADER_DEFINED 1
#include <dl.h>

vtkLibHandle vtkDynamicLoader::OpenLibrary(const char* libname )
{
  return shl_load(libname, BIND_DEFERRED | DYNAMIC_PATH, 0L);
}

int vtkDynamicLoader::CloseLibrary(vtkLibHandle lib)
{
  return 0;
}

void* vtkDynamicLoader::GetSymbolAddress(vtkLibHandle lib, const char* sym)
{ 
  void* addr;
  int status;
  
  status = shl_findsym (&lib, sym, TYPE_PROCEDURE, &addr);
  return (status < 0) ? (void*)0 : addr;
}

const char* vtkDynamicLoader::LibPrefix()
{ 
  return "lib";
}

const char* vtkDynamicLoader::LibExtension()
{
  return ".sl";
}

const char* vtkDynamicLoader::LastError()
{
  return 0;
}
#endif



// ---------------------------------------------------------------
// 2. Implementation for Darwin (including OSX) Machines

#ifdef __APPLE__
#define VTKDYNAMICLOADER_DEFINED
#include <mach-o/dyld.h>

vtkLibHandle vtkDynamicLoader::OpenLibrary(const char* libname )
{
  NSObjectFileImageReturnCode rc;
  NSObjectFileImage image;

  rc = NSCreateObjectFileImageFromFile(libname, &image);
  return NSLinkModule(image, libname, TRUE);
}

int vtkDynamicLoader::CloseLibrary(vtkLibHandle lib)
{
  return 0;
}

void* vtkDynamicLoader::GetSymbolAddress(vtkLibHandle lib, const char* sym)
{
    void *result=0;
    if(NSIsSymbolNameDefined(sym)){
         cout << sym << " is defined!" << endl;
         NSSymbol symbol= NSLookupAndBindSymbol(sym);
         if(symbol){
                result = NSAddressOfSymbol(symbol);
         }
  }else{
        cout << sym << " is not defined!" << endl;
 }
  return result;
}

const char* vtkDynamicLoader::LibPrefix()
{
  return "";
}

const char* vtkDynamicLoader::LibExtension()
{
  return ".dylib";
}

const char* vtkDynamicLoader::LastError()
{
  return 0;
}

#endif




// ---------------------------------------------------------------
// 3. Implementation for Windows win32 code
#ifdef _WIN32
#include <windows.h>
#define VTKDYNAMICLOADER_DEFINED 1

vtkLibHandle vtkDynamicLoader::OpenLibrary(const char* libname )
{
#ifdef UNICODE
        wchar_t *libn = new wchar_t [mbstowcs(NULL, libname, 32000)];
        mbstowcs(libn, libname, 32000);
        vtkLibHandle ret = LoadLibrary(libn);
        delete [] libn;
        return ret;
#else
        return LoadLibrary(libname);
#endif
}

int vtkDynamicLoader::CloseLibrary(vtkLibHandle lib)
{
  return (int)FreeLibrary(lib);
}

void* vtkDynamicLoader::GetSymbolAddress(vtkLibHandle lib, const char* sym)
{ 
#ifdef UNICODE
        wchar_t *wsym = new wchar_t [mbstowcs(NULL, sym, 32000)];
        mbstowcs(wsym, sym, 32000);
        void *ret = GetProcAddress(lib, wsym);
        delete [] wsym;
        return ret;
#else
  return GetProcAddress(lib, sym);
#endif
}

const char* vtkDynamicLoader::LibPrefix()
{ 
  return "";
}

const char* vtkDynamicLoader::LibExtension()
{
  return ".dll";
}

const char* vtkDynamicLoader::LastError()
{
  LPVOID lpMsgBuf;

  FormatMessage( 
                FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM,
                NULL,
                GetLastError(),
                MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT), // Default language
                (LPTSTR) &lpMsgBuf,
                0,
                NULL 
                );
  
  // Free the buffer.
  LocalFree( lpMsgBuf );
  static char* str = 0;
  delete [] str;
  str = strcpy(new char[strlen((char*)lpMsgBuf)+1], (char*)lpMsgBuf);
  return str;
}

#endif

// ---------------------------------------------------------------
// 4. Implementation for default UNIX machines.
// if nothing has been defined then use this
#ifndef VTKDYNAMICLOADER_DEFINED
#define VTKDYNAMICLOADER_DEFINED
// Setup for most unix machines
#include <dlfcn.h>

vtkLibHandle vtkDynamicLoader::OpenLibrary(const char* libname )
{
  return dlopen(libname, RTLD_LAZY);
}

int vtkDynamicLoader::CloseLibrary(vtkLibHandle lib)
{
  return (int)dlclose(lib);
}

void* vtkDynamicLoader::GetSymbolAddress(vtkLibHandle lib, const char* sym)
{ 
  return dlsym(lib, sym);
}

const char* vtkDynamicLoader::LibPrefix()
{ 
  return "lib";
}

const char* vtkDynamicLoader::LibExtension()
{
  return ".so";
}

const char* vtkDynamicLoader::LastError()
{
  return dlerror(); 
}
#endif
