/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkWin32ImageWindow.h
  Language:  C++
  Date:      $Date$
  Version:   $Revision$
  Thanks:    Thanks to Matt Turek who developed this class.

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
// .NAME vtkWin32ImageWindow - 2D display window for Windows
// .SECTION Description
// vtkWin32ImageWindow is a concrete subclass of vtkImageWindow.
// It handles 2D rendering under windows.

// .SECTION See Also
// vtkImageWindow

#ifndef __vtkWin32ImageWindow_h
#define __vtkWin32ImageWindow_h


#include 	"vtkImageWindow.h"

class VTK_EXPORT vtkWin32ImageWindow : public vtkImageWindow 
{
public:
  HINSTANCE ApplicationInstance;
  HPALETTE  Palette;
  HDC       DeviceContext;
  HWND      WindowId;
  HWND      ParentId;

  static vtkWin32ImageWindow *New();
  const char *GetClassName() {return "vtkWin32ImageWindow";};
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Swap the front and back buffers. Normally not called by the user.
  void SwapBuffers();

  // output to the viewer.
  vtkWin32ImageWindow *GetOutput(){return this;};
  
  //BTX

  // Description:
  // Set/Get the window id and parent window id.
  HWND GetWindowId(); 
  void SetWindowId(void* id) {this->WindowId = (HWND) id;};
  void SetParentId(void* id) {this->ParentId = (HWND) id;};
  void SetWindowId(HWND);
  void SetParentId(HWND);

  void SetDeviceContext(void* dc) {this->DeviceContext = (HDC) dc;};
  void SetDeviceContext(HDC);
  void SetDisplayId(void *foo) {vtkDebugMacro(<<"SetDisplayID not implemented");};

  void *GetGenericDisplayId() 
        {vtkDebugMacro(<<"Display ID not implemented in Win32."); return (void*) NULL;};
  void *GetGenericWindowId() {return (void*) this->WindowId;};
  void *GetGenericParentId() {return (void*) this->ParentId;};
  void *GetGenericContext() {return (void*) this->DeviceContext;};
  //ETX

  // Description:
  // Set/Get the current size of the window.
  void   SetSize(int,int);
  int   *GetSize();

  // Description:
  // Set/Get the position in screen coordinates of the window.
  int   *GetPosition();
  void   SetPosition(int,int);

  // Description:
  // Set the desired background color for the window.
  void SetBackgroundColor(float r, float g, float b);

  // Description:
  // Erase the window. Normally nor called by the user.
  void EraseWindow();

  unsigned char *GetDIBPtr();
  unsigned char *GetPixelData(int x1, int y1, int x2, int y2, int);
  
  // Description:
  // Creates a Win32 window or sets up an existing window.
  void MakeDefaultWindow();  

  // Description:
  // These methods can be used by MFC applications 
  // to support print preview and printing, or more
  // general rendering into memory. 
  void SetupMemoryRendering(int x, int y, HDC prn);
  void ResumeScreenRendering();
  HDC GetMemoryDC();
  unsigned char *GetMemoryData(){return this->MemoryData;};

protected:
  vtkWin32ImageWindow();
  ~vtkWin32ImageWindow();
  vtkWin32ImageWindow(const vtkWin32ImageWindow&) {};
  void operator=(const vtkWin32ImageWindow&) {};

  // the following is used to support rendering into memory
  BITMAPINFO MemoryDataHeader;
  HBITMAP MemoryBuffer;
  unsigned char *MemoryData;	// the data in the DIBSection
  HDC MemoryHdc;
  int ScreenMapped;
  int ScreenWindowSize[2];
  HDC ScreenDeviceContext;

  int OwnWindow; // do we create this window ?

  unsigned char *DIBPtr;	// the data in the DIBSection
  int SwapFlag;
  HDC CompatHdc;
  HDC OldHdc;
  HBITMAP BackBuffer;
  BITMAPINFO DataHeader;
};

#endif
