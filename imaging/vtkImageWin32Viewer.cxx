/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkImageWin32Viewer.cxx
  Language:  C++
  Date:      $Date$
  Version:   $Revision$


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
#include "vtkImageWin32Viewer.h"

//----------------------------------------------------------------------------
/* 
 * This templated routine calculates effective lover and upper limits 
 * for a window of values of type T, lower and upper. 
 * It also returns in variable hit the flag with value 1, if the 
 * interval [f_lower, f_upper) is above the type range, -1,  if the 
 * interval [f_lower, f_upper) is below the type range, and 0, if 
 * there is intersection between this interval and the type range.
 */
template <class T>
static void clamps ( vtkImageRegion *region, float w, float l, T& lower, T& upper, int& hit)
{
  float f_lower, f_upper;
  float range[2];

  region->GetData()->GetScalars()->GetDataTypeRange( range );
  
  f_lower = l - fabs(w) / 2.0;
  f_upper = f_lower + fabs(w);

  // Look if we above the type range
  if ( f_lower > range[1] )
    {
      hit = 1;
      return;
    }
  // Look if we below the type range
  if ( f_upper <= range[0] )
    {
      hit = -1;
      return;
    }

  // Set the correct lower value
   if ( f_lower >= range[0] )
    {
      lower = (T) f_lower;
    }
  else
    {
      lower = (T) range[0];
    }

  // Set the correct upper value
  if ( f_upper <= range[1] )
    {
      upper = (T) f_upper;
    }
  else
    {
      upper = (T) range[1];
    }

  hit = 0;
}

vtkImageWin32Viewer::vtkImageWin32Viewer()
{
  this->ApplicationInstance = NULL;
  this->Palette = NULL;
  this->WindowId = 0;
  this->ParentId = 0;
  this->DeviceContext = (HDC)0;
  this->OwnWindow = 0;
  this->HBitmap = (HBITMAP)0;
  this->DataOut = NULL;

  if ( this->WindowName )
    delete [] this->WindowName;
  this->WindowName = strdup("Visualization Toolkit - ImageWin32");
}


//----------------------------------------------------------------------------
vtkImageWin32Viewer::~vtkImageWin32Viewer()
{
  if (this->HBitmap)
    {
    DeleteObject(this->HBitmap);
    this->HBitmap = (HBITMAP)0;
    }
}


//----------------------------------------------------------------------------
void vtkImageWin32Viewer::PrintSelf(ostream& os, vtkIndent indent)
{
  vtkImageViewer::PrintSelf(os, indent);
}

// Description:
// Get the position in screen coordinates of the window.
int *vtkImageWin32Viewer::GetPosition(void)
{
  // if we aren't mapped then just return the ivar 
  if (!this->Mapped)
    {
    return(this->Position);
    }

  //  Find the current window position 
//  x,y,&this->Position[0],&this->Position[1],&child);

  return this->Position;
}

void vtkImageWin32Viewer::SetPosition(int x, int y)
{
  static int resizing = 0;

  if ((this->Position[0] != x) || (this->Position[1] != y))
    {
    this->Modified();
    this->Position[0] = x;
    this->Position[1] = y;
    if (this->Mapped)
      {
      if (!resizing)
        {
        resizing = 1;
   
        SetWindowPos(this->WindowId,HWND_TOP,x,y,
            0, 0, SWP_NOSIZE | SWP_NOZORDER);
        resizing = 0;
        }
      }
    }
}


//----------------------------------------------------------------------------
// Description:
// A templated function that handles gray scale images.
template <class T>
static void vtkImageWin32ViewerRenderGray(vtkImageWin32Viewer *self, 
					  vtkImageRegion *region,
					  T *inPtr, unsigned char *outPtr,
					  float shift, float scale)
{
  unsigned char colorIdx;
  T *inPtr0, *inPtr1, *endPtr;
  int inMin0, inMax0, inMin1, inMax1;
  int inInc0, inInc1;
  int idx1;
  T   lower, upper;
  int rowAdder;
  int *Size;
  int hit;
  int sscale;
  int i_lower;

  clamps ( region, self->GetColorWindow(), self->GetColorLevel(), 
	   lower, upper, hit );
  i_lower = -shift;

  // Selection of the constant 4096.0 may be changed in the future
  // or may be type dependent
  sscale = scale*4096.0;

  region->GetExtent(inMin0, inMax0, inMin1, inMax1);
  region->GetIncrements(inInc0, inInc1);

  // crop to current window size
  Size = self->GetSize();
  if ((inMax0 - inMin0) >= Size[0])
  {
	  inMax0 =  inMin0 - 1 + Size[0];
  }
  if ((inMax1 - inMin1) >= Size[1])
  {
	  inMax1 =  inMin1 - 1 + Size[1];
  }

  // Loop through in regions pixels
  rowAdder = (4 - ((inMax0-inMin0 + 1)*3)%4)%4;
  inPtr1 = inPtr;
  if ( hit == 1 ) 
    {
      // type range is to the left of the window
      for (idx1 = inMin1; idx1 <= inMax1; idx1++)
	{
	  inPtr0 = inPtr1;
	  endPtr = inPtr0 + inInc0*(inMax0 - inMin0 + 1);
	  while (inPtr0 != endPtr)
	    {
	      *outPtr++ = 0;
	      *outPtr++ = 0;
	      *outPtr++ = 0;
	      inPtr0 += inInc0;
	    }
	  // rows must be a multiple of four bytes
	  // so pad it if neccessary
	  outPtr += rowAdder;
	  inPtr1 += inInc1;
	} 
    }
  else if ( hit == -1 ) 
    {
      // type range is to the right of the window
      for (idx1 = inMin1; idx1 <= inMax1; idx1++)
	{
	  inPtr0 = inPtr1;
	  endPtr = inPtr0 + inInc0*(inMax0 - inMin0 + 1);
	  while (inPtr0 != endPtr)
	    {
	      *outPtr++ = 255;
	      *outPtr++ = 255;
	      *outPtr++ = 255;
	      inPtr0 += inInc0;
	    }
	  // rows must be a multiple of four bytes
	  // so pad it if neccessary
	  outPtr += rowAdder;
	  inPtr1 += inInc1;
	} 
    }
  else // type range intersects with the window
    {
      for (idx1 = inMin1; idx1 <= inMax1; idx1++)
	{
	  inPtr0 = inPtr1;
	  endPtr = inPtr0 + inInc0*(inMax0 - inMin0 + 1);
	  while (inPtr0 != endPtr)
	    {
	      if (*inPtr0 < lower) 
		{
		  *outPtr++ = 0;
		  *outPtr++ = 0;
		  *outPtr++ = 0;
		}
	      else if (*inPtr0 >= upper)
		{
		  *outPtr++ = 255;
		  *outPtr++ = 255;
		  *outPtr++ = 255;
		}
	      else
		{
		  // Because i_lower and sscale are of integer type
		  // this is fast for all types used by this
		  // template (float is treated separately).
                  colorIdx = (*inPtr0 - i_lower) * sscale >> 12;
		  *outPtr++ = colorIdx;
		  *outPtr++ = colorIdx;
		  *outPtr++ = colorIdx;
		}
	      inPtr0 += inInc0;
	    }
	  // rows must be a multiple of four bytes
	  // so pad it if neccessary
	  outPtr += rowAdder;
	  inPtr1 += inInc1;
	}
    }
}

//----------------------------------------------------------------------------
// Description:
// A function that handles gray scale images with float data.
static void vtkImageWin32ViewerRenderFloatGray(vtkImageWin32Viewer *self, 
					  vtkImageRegion *region,
					  float *inPtr, unsigned char *outPtr,
					  float shift, float scale)
{
  unsigned char colorIdx;
  float *inPtr0, *inPtr1, *endPtr;
  int inMin0, inMax0, inMin1, inMax1;
  int inInc0, inInc1;
  int idx1;
  float  lower, upper;
  int rowAdder;
  int *Size;

  lower = -shift;
  upper = lower + 255.0/scale;

  region->GetExtent(inMin0, inMax0, inMin1, inMax1);
  region->GetIncrements(inInc0, inInc1);

  // crop to current window size
  Size = self->GetSize();
  if ((inMax0 - inMin0) >= Size[0])
  {
	  inMax0 =  inMin0 - 1 + Size[0];
  }
  if ((inMax1 - inMin1) >= Size[1])
  {
	  inMax1 =  inMin1 - 1 + Size[1];
  }

  // Loop through in regions pixels
  rowAdder = (4 - ((inMax0-inMin0 + 1)*3)%4)%4;
  inPtr1 = inPtr;
  for (idx1 = inMin1; idx1 <= inMax1; idx1++)
    {
      inPtr0 = inPtr1;
      endPtr = inPtr0 + inInc0*(inMax0 - inMin0 + 1);
      while (inPtr0 != endPtr)
	{
	  if (*inPtr0 < lower) 
	    {
	      *outPtr++ = 0;
	      *outPtr++ = 0;
	      *outPtr++ = 0;
	    }
	  else if (*inPtr0 >= upper)
	    {
	      *outPtr++ = 255;
	      *outPtr++ = 255;
	      *outPtr++ = 255;
	    }
	  else
	    {
	      colorIdx = (*inPtr0 - lower) * scale ;
	      *outPtr++ = colorIdx;
	      *outPtr++ = colorIdx;
	      *outPtr++ = colorIdx;
	    }
	  inPtr0 += inInc0;
	}
      // rows must be a multiple of four bytes
      // so pad it if neccessary
      outPtr += rowAdder;
      inPtr1 += inInc1;
    }
}

//----------------------------------------------------------------------------
// Description:
// A templated function that handles color images. (only True Color 24 bit)
template <class T>
static void vtkImageWin32ViewerRenderColor(vtkImageWin32Viewer *self, 
					   vtkImageRegion *region,
					   T *redPtr, T *greenPtr, T *bluePtr,
					   unsigned char *outPtr,
					   float shift, float scale)
{
  unsigned char red, green, blue;
  T *redPtr0, *redPtr1;
  T *bluePtr0, *bluePtr1;
  T *greenPtr0, *greenPtr1;
  int inMin0, inMax0, inMin1, inMax1;
  int inInc0, inInc1;
  int idx0, idx1;
  T   lower, upper;
  int rowAdder;
  
  region->GetExtent(inMin0, inMax0, inMin1, inMax1);
  region->GetIncrements(inInc0, inInc1);
  
  lower = -shift;
  upper = lower + 255.0/scale;

  // Loop through in regions pixels
  redPtr1 = redPtr;
  greenPtr1 = greenPtr;
  bluePtr1 = bluePtr;
  rowAdder = (4 - ((inMax0-inMin0 + 1)*3)%4)%4;
  for (idx1 = inMin1; idx1 <= inMax1; idx1++)
    {
    redPtr0 = redPtr1;
    greenPtr0 = greenPtr1;
    bluePtr0 = bluePtr1;
    for (idx0 = inMin0; idx0 <= inMax0; idx0++)
      {
      if (*redPtr0 <= lower) red = 0;
      else if (*redPtr0 >= upper) red = 255;
      else red = (unsigned char)(((float)(*redPtr0) + shift) * scale);

      if (*greenPtr0 <= lower) green = 0;
      else if (*greenPtr0 >= upper) green = 255;
      else green = (unsigned char)(((float)(*greenPtr0) + shift) * scale);
  
      if (*bluePtr0 <= lower) blue = 0;
      else if (*bluePtr0 >= upper) blue = 255;
      else blue = (unsigned char)(((float)(*bluePtr0) + shift) * scale);
      *outPtr++ = blue;
      *outPtr++ = green;
      *outPtr++ = red;

      redPtr0 += inInc0;
      greenPtr0 += inInc0;
      bluePtr0 += inInc0;
      }
    // rows must be a multiple of four bytes
    // so pad it if neccessary
    outPtr += rowAdder;

    redPtr1 += inInc1;
    greenPtr1 += inInc1;
    bluePtr1 += inInc1;
    }
}

void vtkImageWin32Viewer::SetSize(int x, int y)
{
  static int resizing = 0;

  if ((this->Size[0] != x) || (this->Size[1] != y))
    {
    this->Modified();
    this->Size[0] = x;
    this->Size[1] = y;
    if (this->Mapped)
      {
      if (!resizing)
        {
        BITMAPINFO dataHeader;
        resizing = 1;

        if (this->ParentId)
          {
          SetWindowPos(this->WindowId,HWND_TOP,0,0,
            x, y, SWP_NOMOVE | SWP_NOZORDER);
          }
        else
          {
          SetWindowPos(this->WindowId,HWND_TOP,0,0,
            x+2*GetSystemMetrics(SM_CXFRAME),
            y+2*GetSystemMetrics(SM_CYFRAME) +GetSystemMetrics(SM_CYCAPTION),
            SWP_NOMOVE | SWP_NOZORDER);
          }

        // free and then realloc the DIB
        if (this->HBitmap)
          {
          DeleteObject(this->HBitmap);
          }
        dataHeader.bmiHeader.biSize = 40;
        dataHeader.bmiHeader.biWidth = x;
        dataHeader.bmiHeader.biHeight = y;
        dataHeader.bmiHeader.biPlanes = 1;
        dataHeader.bmiHeader.biBitCount = 24;
        dataHeader.bmiHeader.biCompression = BI_RGB;
        dataHeader.bmiHeader.biSizeImage = ((x*3+3)/4)*4*y;
        dataHeader.bmiHeader.biClrUsed = 0;
        dataHeader.bmiHeader.biClrImportant = 0;  

        // try using a DIBsection
        this->HBitmap = CreateDIBSection(this->DeviceContext, &dataHeader, 
                             DIB_RGB_COLORS, (void **)(&(this->DataOut)), NULL, 0);
        resizing = 0;
        }
      }
    }
}

// Description:
// Get the current size of the window.
int *vtkImageWin32Viewer::GetSize(void)
{
  RECT rect;

  // if we aren't mapped then just return the ivar 
  if (!this->Mapped)
    {
    return(this->Size);
    }

  //  Find the current window size 
  GetClientRect(this->WindowId, &rect);

  this->Size[0] = rect.right;
  this->Size[1] = rect.bottom;
  
  return this->Size;
}

//----------------------------------------------------------------------------
// Expects region to be X, Y, components
void vtkImageWin32Viewer::RenderRegion(vtkImageRegion *region)
{
  int dataWidth, width, height;
  int size;
  int extent[6];
  unsigned char *dataOut;
  void *ptr0, *ptr1, *ptr2;
  float shift, scale;
  HDC compatDC;
  HBITMAP hOldBitmap;

  if ( ! region)
    {
    // open the window anyhow if one has not been set.
    if (!this->DeviceContext)
      {
      // use default size if not specified
      if (this->Size[0] == 0 )
        {
        this->Size[0] = 256;
        this->Size[1] = 256;
        }
      this->MakeDefaultWindow();
      }
    return;
    }
  
  // Determine the size of the displayed region.
  region->GetExtent(3, extent);
  width = (extent[1] - extent[0] + 1);
  height = (extent[3] - extent[2] + 1);

  // In case a window has not been set.
  if (!this->DeviceContext)
    {
    // use default size if not specified
    if (this->Size[0] == 0 )
      {
      this->Size[0] = width;
      this->Size[1] = height;
      }
    this->MakeDefaultWindow();
    }
  
  dataWidth = ((width*3+3)/4)*4;
  shift = this->ColorWindow / 2.0 - this->ColorLevel;
  scale = 255.0 / this->ColorWindow;


  BITMAP bm;
  GetObject(this->HBitmap, sizeof (BITMAP), (LPSTR) &bm);

  // vtkDebugMacro(<< "vtkImageWin32Viewer::RenderRegion - Bitmap width: " << bm.bmWidth);
  // vtkDebugMacro(<< "vtkImageWin32Viewer::RenderRegion - Bitmap height: " << bm.bmHeight);

  // create the DIBSection if not done already
  if (!this->HBitmap)
    {
    BITMAPINFO dataHeader;
    dataHeader.bmiHeader.biSize = 40;
    dataHeader.bmiHeader.biWidth = width;
    dataHeader.bmiHeader.biHeight = height;
    dataHeader.bmiHeader.biPlanes = 1;
    dataHeader.bmiHeader.biBitCount = 24;
    dataHeader.bmiHeader.biCompression = BI_RGB;
    dataHeader.bmiHeader.biSizeImage = dataWidth*height;
    dataHeader.bmiHeader.biClrUsed = 0;
    dataHeader.bmiHeader.biClrImportant = 0;  

    // try using a DIBsection
    this->HBitmap = CreateDIBSection(this->DeviceContext, &dataHeader, 
                             DIB_RGB_COLORS, (void **)(&(this->DataOut)), NULL, 0);
    }
   // free and then realloc the DIB if it needs to change size (Window resize)
   // if region size differs from bitmap size, reallocate the bitmap
   else if ((width != bm.bmWidth) || (height != bm.bmHeight))
     {

	// vtkDebugMacro(<< "vtkImageWin32Viewer::RenderRegion - Changing bitmap size to: "
	//   << width << "," << height << "(" << dataWidth*height << " bytes)");
    	
    DeleteObject(this->HBitmap);

    BITMAPINFO dataHeader;
    dataHeader.bmiHeader.biSize = 40;
    dataHeader.bmiHeader.biWidth = width;
    dataHeader.bmiHeader.biHeight = height;
    dataHeader.bmiHeader.biPlanes = 1;
    dataHeader.bmiHeader.biBitCount = 24;
    dataHeader.bmiHeader.biCompression = BI_RGB;
    dataHeader.bmiHeader.biSizeImage = dataWidth*height;
    dataHeader.bmiHeader.biClrUsed = 0;
    dataHeader.bmiHeader.biClrImportant = 0;  

    // try using a DIBsection
    this->HBitmap = CreateDIBSection(this->DeviceContext, &dataHeader, 
                         DIB_RGB_COLORS, (void **)(&(this->DataOut)), NULL, 0);

     }
   
  int min = 0;
  int max = 0;
  int dim = 0;
  region->GetAxisExtent(VTK_IMAGE_COMPONENT_AXIS, min, max);
  dim = max - min + 1;

  if (dim > 1)    
    {

      ptr0 = region->GetScalarPointer(extent[0], extent[2], 
				      this->RedComponent);
      ptr1 = region->GetScalarPointer(extent[0], extent[2], 
				      this->GreenComponent);
      ptr2 = region->GetScalarPointer(extent[0], extent[2], 
				      this->BlueComponent);


    if ( ! ptr0 ||! ptr1 || ! ptr2)
      {
      vtkErrorMacro("Render: Could not get date. Check that RGB are in range");
      return;
      }
  

    // Call the appropriate templated function
    switch (region->GetScalarType())
      {
      case VTK_FLOAT:
	vtkImageWin32ViewerRenderColor(this, region, 
			   (float *)(ptr0),(float *)(ptr1),(float *)(ptr2), 
			   this->DataOut, shift, scale);
	break;
      case VTK_INT:
	vtkImageWin32ViewerRenderColor(this, region, 
			   (int *)(ptr0), (int *)(ptr1), (int *)(ptr2), 
			   this->DataOut, shift, scale);
	break;
      case VTK_SHORT:
	vtkImageWin32ViewerRenderColor(this, region, 
			   (short *)(ptr0),(short *)(ptr1),(short *)(ptr2), 
			   this->DataOut, shift, scale);
	break;
      case VTK_UNSIGNED_SHORT:
	vtkImageWin32ViewerRenderColor(this, region, (unsigned short *)(ptr0),
			   (unsigned short *)(ptr1),(unsigned short *)(ptr2), 
			    this->DataOut, shift, scale);
	break;
      case VTK_UNSIGNED_CHAR:
	vtkImageWin32ViewerRenderColor(this, region, (unsigned char *)(ptr0), 
			   (unsigned char *)(ptr1),(unsigned char *)(ptr2), 
			    this->DataOut, shift, scale);
	break;
      }
    }
  else
    {
    // GrayScale images.
    ptr0 = region->GetScalarPointer();
    // Call the appropriate templated function
    switch (region->GetScalarType())
      {
      case VTK_FLOAT:
	      vtkImageWin32ViewerRenderFloatGray(this, region, (float *)(ptr0), this->DataOut, shift, scale);
	      break;
      case VTK_INT:
	      vtkImageWin32ViewerRenderGray(this, region, (int *)(ptr0), this->DataOut,
				                              shift, scale);
	      break;
      case VTK_SHORT:
	      vtkImageWin32ViewerRenderGray(this, region, (short *)(ptr0), this->DataOut,
				                                   shift, scale);
	      break;
      case VTK_UNSIGNED_SHORT:
	      vtkImageWin32ViewerRenderGray(this, region, (unsigned short *)(ptr0), 
				                                   this->DataOut, shift, scale);
	      break;
      case VTK_UNSIGNED_CHAR:
	      vtkImageWin32ViewerRenderGray(this, region, (unsigned char *)(ptr0), 
				      this->DataOut, shift, scale);
	      break;
      }   
    }
  
  compatDC = CreateCompatibleDC(this->DeviceContext);  
  hOldBitmap = (HBITMAP)SelectObject(compatDC,this->HBitmap);

#if 0
  // #####

  vtkDebugMacro(<<"vtkWin32ImageMapper::RenderRegion - COLORES: " <<
	              GetDeviceCaps(this->DeviceContext, COLORRES));

  vtkDebugMacro(<<"vtkWin32ImageMapper::RenderRegion - NUMCOLORS: " <<
	              GetDeviceCaps(this->DeviceContext, NUMCOLORS));

  vtkDebugMacro(<<"vtkWin32ImageMapper::RenderRegion - PLANES: " <<
	              GetDeviceCaps(this->DeviceContext, PLANES));

#endif

  BitBlt(this->DeviceContext,0,0,width,height,compatDC,0,0,SRCCOPY);
  SelectObject(compatDC, hOldBitmap);
  DeleteDC(compatDC);
    
}

void vtkImageWin32ViewerSetupRGBPixelFormat(HDC hDC)
{
    PIXELFORMATDESCRIPTOR pfd = {
        sizeof(PIXELFORMATDESCRIPTOR),  /* size */
        1,                              /* version */
        PFD_DRAW_TO_WINDOW,
        PFD_TYPE_RGBA,                   /* color type */
        24,                             /* prefered color depth */
        0, 0, 0, 0, 0, 0,               /* color bits (ignored) */
        0,                              /* no alpha buffer */
        0,                              /* alpha bits (ignored) */
        0,                              /* no accumulation buffer */
        0, 0, 0, 0,                     /* accum bits (ignored) */
        0,                              /* depth buffer */
        0,                              /* no stencil buffer */
        0,                              /* no auxiliary buffers */
        PFD_MAIN_PLANE,                 /* main layer */
        0,                              /* reserved */
        0, 0, 0,                        /* no layer, visible, damage masks */
    };
    int pixelFormat;

    pixelFormat = ChoosePixelFormat(hDC, &pfd);
    if (pixelFormat == 0) {
        MessageBox(WindowFromDC(hDC), "ChoosePixelFormat failed.", "Error",
                MB_ICONERROR | MB_OK);
        exit(1);
    }

    if (SetPixelFormat(hDC, pixelFormat, &pfd) != TRUE) {
        MessageBox(WindowFromDC(hDC), "SetPixelFormat failed.", "Error",
                MB_ICONERROR | MB_OK);
        exit(1);
    }
}

void vtkImageWin32ViewerSetupGrayPixelFormat(HDC hDC)
{
    PIXELFORMATDESCRIPTOR pfd = {
        sizeof(PIXELFORMATDESCRIPTOR),  /* size */
        1,                              /* version */
        PFD_DRAW_TO_WINDOW,
        PFD_TYPE_COLORINDEX,            /* color type */
        8,                              /* prefered color depth */
        0, 0, 0, 0, 0, 0,               /* color bits (ignored) */
        0,                              /* no alpha buffer */
        0,                              /* alpha bits (ignored) */
        0,                              /* no accumulation buffer */
        0, 0, 0, 0,                     /* accum bits (ignored) */
        0,                              /* depth buffer */
        0,                              /* no stencil buffer */
        0,                              /* no auxiliary buffers */
        PFD_MAIN_PLANE,                 /* main layer */
        0,                              /* reserved */
        0, 0, 0,                        /* no layer, visible, damage masks */
    };
    int pixelFormat;

    pixelFormat = ChoosePixelFormat(hDC, &pfd);
    if (pixelFormat == 0) {
        MessageBox(WindowFromDC(hDC), "ChoosePixelFormat failed.", "Error",
                MB_ICONERROR | MB_OK);
        exit(1);
    }

    if (SetPixelFormat(hDC, pixelFormat, &pfd) != TRUE) {
        MessageBox(WindowFromDC(hDC), "SetPixelFormat failed.", "Error",
                MB_ICONERROR | MB_OK);
        exit(1);
    }
}

// creates and applies a RGB palette
void vtkImageWin32ViewerSetupRGBPalette(HDC hDC, 
				     vtkImageWin32Viewer *me)
{
  int pixelFormat = GetPixelFormat(hDC);
  PIXELFORMATDESCRIPTOR pfd;
  LOGPALETTE* pPal;
  int paletteSize;
    
  DescribePixelFormat(hDC, pixelFormat, sizeof(PIXELFORMATDESCRIPTOR), &pfd);
  
  if (pfd.dwFlags & PFD_NEED_PALETTE) 
    {
    paletteSize = 1 << pfd.cColorBits;
    } 
  else 
    {
    return;
    }
  
  pPal = (LOGPALETTE*)
    malloc(sizeof(LOGPALETTE) + paletteSize * sizeof(PALETTEENTRY));
  pPal->palVersion = 0x300;
  pPal->palNumEntries = paletteSize;
  
  /* build a simple RGB color palette */
  {
  int redMask = (1 << pfd.cRedBits) - 1;
  int greenMask = (1 << pfd.cGreenBits) - 1;
  int blueMask = (1 << pfd.cBlueBits) - 1;
  int i;
  
  for (i=0; i<paletteSize; ++i) 
    {
    pPal->palPalEntry[i].peRed =
      (((i >> pfd.cRedShift) & redMask) * 255) / redMask;
    pPal->palPalEntry[i].peGreen =
      (((i >> pfd.cGreenShift) & greenMask) * 255) / greenMask;
    pPal->palPalEntry[i].peBlue =
      (((i >> pfd.cBlueShift) & blueMask) * 255) / blueMask;
    pPal->palPalEntry[i].peFlags = 0;
    }
  }
  
  me->Palette = CreatePalette(pPal);
  free(pPal);
  
  if (me->Palette) 
    {
    SelectPalette(hDC, me->Palette, FALSE);
    RealizePalette(hDC);
    }

}

void vtkImageWin32ViewerSetupGrayPalette(HDC hDC, 
				     vtkImageWin32Viewer *me)
{
  int pixelFormat = GetPixelFormat(hDC);
  PIXELFORMATDESCRIPTOR pfd;
  LOGPALETTE* pPal;
  int paletteSize;
  
  DescribePixelFormat(hDC, pixelFormat, sizeof(PIXELFORMATDESCRIPTOR), &pfd);
  
  // we always want a palette on 8 bit displays
  if (pfd.cColorBits == 8 || pfd.dwFlags & PFD_NEED_PALETTE) 
    {
    paletteSize = 1 << pfd.cColorBits;
    } 
  else 
    {
    return;
    }
  
  pPal = (LOGPALETTE*)
    malloc(sizeof(LOGPALETTE) + paletteSize * sizeof(PALETTEENTRY));
  pPal->palVersion = 0x300;
  pPal->palNumEntries = paletteSize;
  
  /* build a simple RGB color palette */
  {
  int redMask = (1 << pfd.cRedBits) - 1;
  int greenMask = (1 << pfd.cGreenBits) - 1;
  int blueMask = (1 << pfd.cBlueBits) - 1;
  int i;
  
  for (i=0; i<paletteSize; ++i) 
    {
    pPal->palPalEntry[i].peRed = (255*i)/paletteSize;
    pPal->palPalEntry[i].peGreen = (255*i)/paletteSize;
    pPal->palPalEntry[i].peBlue = (255*i)/paletteSize;
    pPal->palPalEntry[i].peFlags = 0;
    }
  }
  
  me->Palette = CreatePalette(pPal);
  free(pPal);
  
  if (me->Palette) 
    {
    SelectPalette(hDC, me->Palette, FALSE);
    RealizePalette(hDC);
    }
}

// used to pass info into the create routine because there doesn't
// seem to be another way. Could be a problem for multithreaded
// apps but this is unlikely since this doesn't get called very
// often at all.
vtkImageWin32Viewer *vtkImageWin32ViewerPtr = NULL;

LRESULT APIENTRY vtkImageWin32ViewerWndProc(HWND hWnd, UINT message, 
					    WPARAM wParam, LPARAM lParam)
{
  vtkImageWin32Viewer *me =   
    (vtkImageWin32Viewer *)GetWindowLong(hWnd,GWL_USERDATA);

  switch (message) 
    {
    case WM_CREATE:
      {
        me = vtkImageWin32ViewerPtr;
        SetWindowLong(hWnd,GWL_USERDATA,(LONG)me);
        me->DeviceContext = GetDC(hWnd);
        if (me->GetGrayScaleHint())
          {
          vtkImageWin32ViewerSetupGrayPixelFormat(me->DeviceContext);
          vtkImageWin32ViewerSetupGrayPalette(me->DeviceContext,me);
          }
        else
          {
          vtkImageWin32ViewerSetupRGBPixelFormat(me->DeviceContext);
          vtkImageWin32ViewerSetupRGBPalette(me->DeviceContext,me);
          }
        return 0;
      }
    case WM_DESTROY:
        if (me->Palette)
          {
          DeleteObject(me->Palette);
          me->Palette = NULL;
          }
        ReleaseDC(me->WindowId, me->DeviceContext);
        return 0;
    case WM_SIZE:
        /* track window size changes */
        if (me->DeviceContext) 
          {
          me->SetSize((int) LOWORD(lParam),(int) HIWORD(lParam));
          return 0;
          }
    case WM_PALETTECHANGED:
        /* realize palette if this is *not* the current window */
        if (me->DeviceContext && me->Palette && (HWND) wParam != hWnd) 
          {
          UnrealizeObject(me->Palette);
          SelectPalette(me->DeviceContext, me->Palette, FALSE);
          RealizePalette(me->DeviceContext);
          me->Render();
          break;
          }
        break;
    case WM_QUERYNEWPALETTE:
        /* realize palette if this is the current window */
        if (me->DeviceContext && me->Palette) 
          {
          UnrealizeObject(me->Palette);
          SelectPalette(me->DeviceContext, me->Palette, FALSE);
          RealizePalette(me->DeviceContext);
          me->Render();
          return TRUE;
          }
        break;
    case WM_PAINT:
        {
        PAINTSTRUCT ps;
        BeginPaint(hWnd, &ps);
        if (me->DeviceContext) 
          {
          me->Render();
          }
        EndPaint(hWnd, &ps);
        return 0;
        }
        break;
    default:
        break;
    }
    return DefWindowProc(hWnd, message, wParam, lParam);
}


//----------------------------------------------------------------------------
void vtkImageWin32Viewer::MakeDefaultWindow() 
{
  static int count = 0;

  // create our own window if not already set
  this->OwnWindow = 0;

  // get the applicaiton instance if we don't have one already
  if (!this->ApplicationInstance)
    {
    // if we have a parent window get the app instance from it
    if (this->ParentId)
      {
      this->ApplicationInstance = 
	(HINSTANCE)GetWindowLong(this->ParentId,GWL_HINSTANCE);
      }
    else
      {
      this->ApplicationInstance = GetModuleHandle(NULL);
      }
    }
  if (!this->WindowId)
    {
    WNDCLASS wndClass;
    
    if(this->WindowName) delete [] this->WindowName;
    int len = strlen( "Visualization Toolkit - ImageWin32 #") 
      + (int)ceil( (double) log10( (double)(count+1) ) ) + 1; 
    this->WindowName = new char [ len ];
    sprintf(this->WindowName,"Visualization Toolkit - ImageWin32 #%i",count++);
    
    // has the class been registered ?
    if (!GetClassInfo(this->ApplicationInstance,"vtkImage",&wndClass))
        {
        wndClass.style = CS_HREDRAW | CS_VREDRAW;
        wndClass.lpfnWndProc = vtkImageWin32ViewerWndProc;
        wndClass.cbClsExtra = 0;
        wndClass.cbWndExtra = 0;
        wndClass.hInstance = this->ApplicationInstance;
        wndClass.hIcon = LoadIcon(NULL, IDI_APPLICATION);
        wndClass.hCursor = LoadCursor(NULL, IDC_ARROW);
        wndClass.hbrBackground = (HBRUSH)GetStockObject(BLACK_BRUSH);
        wndClass.lpszMenuName = NULL;
        wndClass.lpszClassName = "vtkImage";
        RegisterClass(&wndClass);
        }
    
    /* create window */
    // use poor mans mutex
    if (vtkImageWin32ViewerPtr)
      {
      vtkErrorMacro("Two windows being created at the same time");
      }
    vtkImageWin32ViewerPtr = this;
    if (this->ParentId)
      {
      this->WindowId = 
	      CreateWindow("vtkImage", this->WindowName,
		     WS_CHILD | WS_CLIPCHILDREN | WS_CLIPSIBLINGS,
		     0, 0, this->Size[0], this->Size[1],
		     this->ParentId, NULL, this->ApplicationInstance, NULL);
      }
    else
      {
      this->WindowId = 
	      CreateWindow("vtkImage", this->WindowName,
		     WS_OVERLAPPEDWINDOW | WS_CLIPCHILDREN | WS_CLIPSIBLINGS,
		     0, 0, this->Size[0]+2*GetSystemMetrics(SM_CXFRAME),
         this->Size[1] + 2*GetSystemMetrics(SM_CYFRAME) + 
         GetSystemMetrics(SM_CYCAPTION),
		     NULL, NULL, this->ApplicationInstance, NULL);
      }
    vtkImageWin32ViewerPtr = NULL;
    if (!this->WindowId)
      {
      vtkErrorMacro("Could not create window, error:  " << GetLastError());
      return;
      }
    
    /* display window */
    ShowWindow(this->WindowId, SW_SHOW);
    this->OwnWindow = 1;
    }
  this->Mapped = 1;
}

// Description:
// Get the window id.
HWND vtkImageWin32Viewer::GetWindowId()
{
  vtkDebugMacro(<< "Returning WindowId of " << this->WindowId << "\n"); 

  return this->WindowId;
}

// Description:
// Set the window id to a pre-existing window.
void vtkImageWin32Viewer::SetWindowId(HWND arg)
{
  vtkDebugMacro(<< "Setting WindowId to " << arg << "\n"); 

  this->WindowId = arg;
}

// Description:
// Set the window id to a pre-existing window.
void vtkImageWin32Viewer::SetParentId(HWND arg)
{
  vtkDebugMacro(<< "Setting ParentId to " << arg << "\n"); 

  this->ParentId = arg;
}

