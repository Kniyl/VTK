#include "vtkXImageMapper.h"
#include "vtkImageWindow.h"
#include "vtkXImageWindow.h"
#include "vtkProperty2D.h"

vtkXImageMapper::vtkXImageMapper()
{
  this->Image = (XImage*) NULL;
  this->DataOut = NULL;
  this->DataOutSize = 0;
  this->NumberOfColors = 150;
}

vtkXImageMapper::~vtkXImageMapper()
{
  // clean up the render buffer
  if (this->DataOut) delete[] this->DataOut;
  
}

//----------------------------------------------------------------------------
/* 
 * This templated routine calculates effective lower and upper limits 
 * for a window of values of type T, lower and upper. 
 * It also returns in variable hit the flag with value 1, if the 
 * interval [f_lower, f_upper) is above the type range, -1,  if the 
 * interval [f_lower, f_upper) is below the type range, and 0, if 
 * there is intersection between this interval and the type range.
 */
template <class T>
static void vtkXImageMapperClamps ( vtkImageData *data, float w, 
				    float l, T& lower, T& upper, 
				    unsigned char &lower_val, 
				    unsigned char &upper_val)
{
  double f_lower, f_upper, f_lower_val, f_upper_val;
  double range[2];

  data->GetPointData()->GetScalars()->GetDataTypeRange( range );

  f_lower = l - fabs(w) / 2.0;
  f_upper = f_lower + fabs(w);

  // Set the correct lower value
  if ( f_lower <= range[1])
    {
    if (f_lower >= range[0])
      {
      lower = (T) f_lower;
      }
    else
      {
      lower = (T) range[0];
      }
    }
  else
    {
    lower = (T) range[1];
    }
  
  
  // Set the correct upper value
  if ( f_upper >= range[0])
    {
    if (f_upper <= range[1])
      {
      upper = (T) f_upper;
      }
    else
      {
      upper = (T) range[1];
      }
    }
  else
    {
    upper = (T) range [0];
    }
  
  // now compute the lower and upper values
  if (w >= 0)
    {
    f_lower_val = 255.0*(lower - f_lower)/w;
    f_upper_val = 255.0*(upper - f_lower)/w;
    }
  else
    {
    f_lower_val = 255.0 + 255.0*(lower - f_lower)/w;
    f_upper_val = 255.0 + 255.0*(upper - f_lower)/w;
    }
  
  if (f_upper_val > 255) 
    {
    upper_val = 255;
    }
  else if (f_upper_val < 0)
    {
    upper_val = 0;
    }
  else
    {
    upper_val = (unsigned char)(f_upper_val);
    }
  
  if (f_lower_val > 255) 
    {
    lower_val = 255;
    }
  else if (f_lower_val < 0)
    {
    lower_val = 0;
    }
  else
    {
    lower_val = (unsigned char)(f_lower_val);
    }  
}

int vtkXImageMapper::GetCompositingMode(vtkActor2D* actor)
{

  vtkProperty2D* tempProp = actor->GetProperty();
  int compositeMode = tempProp->GetCompositingOperator();

  switch (compositeMode)
  {
  case VTK_BLACK:
	  return GXclear;
  case VTK_NOT_DEST:
	  return GXinvert;
  case VTK_SRC_AND_DEST:
	  return GXand;
  case VTK_SRC_OR_DEST:
	  return  GXor;
  case VTK_NOT_SRC:
	  return GXcopyInverted;
  case VTK_SRC_XOR_DEST:
	  return GXxor;
  case VTK_SRC_AND_notDEST:
	  return GXandReverse;
  case VTK_SRC:
	  return GXcopy;
  case VTK_WHITE:
	  return GXset;
  default:
	  return GXcopy;
  }

}

int vtkXImageMapper::GetXWindowDepth(vtkWindow* window)
{
  Window windowID = (Window) window->GetGenericWindowId();
  Display* displayID = (Display*) window->GetGenericDisplayId();

  XWindowAttributes winAttribs;
  XGetWindowAttributes(displayID, windowID, &winAttribs);

  vtkDebugMacro(<<"vtkXImageMapper::GetXWindowDepth - Returning window depth of: " << winAttribs.depth);
  return winAttribs.depth;

}

void vtkXImageMapper::GetXWindowVisualId(vtkWindow* window, Visual* visualID)
{

  XWindowAttributes winAttribs;

  Window windowID = (Window) window->GetGenericWindowId();
  Display* displayID = (Display*) window->GetGenericDisplayId();
  XGetWindowAttributes(displayID, windowID, &winAttribs);
  *visualID = *winAttribs.visual;

}

int vtkXImageMapper::GetXWindowVisualClass(vtkWindow* window)
{

  Window windowID = (Window) window->GetGenericWindowId();
  Display* displayID = (Display*) window->GetGenericDisplayId();

  // Get the visual class
  int nvisuals = 0;
  XVisualInfo templ;  

  XWindowAttributes winAttribs;
  XGetWindowAttributes(displayID, windowID, &winAttribs);

  templ.visualid = winAttribs.visual->visualid;
  XVisualInfo* visuals = XGetVisualInfo(displayID,
                           VisualIDMask,
                           &templ, &nvisuals);
  if (nvisuals == 0)
    {
    vtkErrorMacro(<< "Could not get visual class");
    }

  vtkDebugMacro(<< "Got visual class");

  int visClass = visuals->c_class;
  XFree(visuals);
  return visClass;

  // return visuals->c_class;
}

void vtkXImageMapper::GetXWindowColorMasks(vtkWindow *window, unsigned long *rmask,
					   unsigned long *gmask, unsigned long *bmask)
{
  Window windowID = (Window) window->GetGenericWindowId();
  Display* displayID = (Display*) window->GetGenericDisplayId();

  XWindowAttributes winAttribs;
  XGetWindowAttributes(displayID, windowID, &winAttribs);
 
  XVisualInfo temp1;
  temp1.visualid = winAttribs.visual->visualid;

  int nvisuals = 0;
  XVisualInfo* visuals = XGetVisualInfo(displayID, VisualIDMask, &temp1, &nvisuals);   

  if (nvisuals == 0)  vtkErrorMacro(<<"Could not get color masks");

  *rmask = visuals->red_mask;
  *gmask = visuals->green_mask;
  *bmask = visuals->blue_mask;

  XFree(visuals);

  return;

}
void vtkXImageMapper::GetXWindowColors(vtkWindow* window, XColor colors[], int ncolors)
{
  int idx = 0;

  Window windowID = (Window) window->GetGenericWindowId();
  Display* displayID = (Display*) window->GetGenericDisplayId();

  XWindowAttributes winAttribs;
  XGetWindowAttributes(displayID, windowID, &winAttribs);
  Colormap colorMap = winAttribs.colormap;

  // Get the colors in the current color map.
  for (idx = 0 ; idx < ncolors; idx++) 
    {
    colors[idx].pixel = idx;
    }

  XQueryColors(displayID, colorMap, colors, ncolors);
  vtkDebugMacro(<<"Got colors.");

}



//----------------------------------------------------------------------------
// Description:
// A templated function that handles gray scale images.
template <class T>
static void vtkXImageMapperRenderGray(vtkXImageMapper *mapper, 
				      vtkViewport *viewport,
                                      vtkImageData *data,
                                      T *inPtr, unsigned char *outPtr)
{

  T *inPtr0, *inPtr1, *endPtr;
  int inMin0, inMax0, inMin1, inMax1;
  int inInc0, inInc1;
  int idx1;
  float shift, scale;
  int colorsMax;
  int visualDepth, visualClass;
  T lower, upper;
  unsigned char lowerPixel, upperPixel, temp;
  int colorIdx;
  unsigned char lower_val, upper_val;
  
  vtkWindow*  window = viewport->GetVTKWindow();

  visualClass = mapper->GetXWindowVisualClass(window);
  visualDepth = mapper->GetXWindowDepth(window);

  XColor *colors = new XColor[256];

  mapper->GetXWindowColors(window, colors, 256);

  shift = mapper->GetColorShift();
  scale = mapper->GetColorScale();

  int* tempExt = mapper->GetInput()->GetUpdateExtent();
  inMin0 = tempExt[0];
  inMax0 = tempExt[1];
  inMin1 = tempExt[2];
  inMax1 = tempExt[3];

  int* tempIncs = data->GetIncrements();
  inInc0 = tempIncs[0];
  inInc1 = tempIncs[1];

  unsigned long rmask = 0;
  unsigned long gmask = 0;
  unsigned long bmask = 0;

  mapper->GetXWindowColorMasks(window, &rmask, &gmask, &bmask);

  // Compute the shifts needed to align the color masks with the 
  // pixels
  
  int rshift = 0;

  while ( ((rmask & 0x80000000) == 0) && (rshift < 32))
    {
    rmask = rmask << 1;
    rshift++;
    }

  int gshift = 0;

  while ( ((gmask & 0x80000000) == 0) && (gshift < 32))
    {
    gmask = gmask << 1;
    gshift++;
    }

  int bshift = 0;

  while ( ((bmask & 0x80000000) == 0) && (bshift < 32))
    {
    bmask = bmask << 1;
    bshift++;
    }

  unsigned long* ulOutPtr = (unsigned long*) outPtr;

  vtkXImageMapperClamps( data, mapper->GetColorWindow(),
			 mapper->GetColorLevel(), 
			 lower, upper, lower_val, upper_val);

  // compute colorsMax, lower and upper pixels.
  if (visualClass == TrueColor)
    {
    colorsMax = upper_val; 
    upperPixel = colorsMax;
    lowerPixel = lower_val;
    }
  else 
    {
    colorsMax = mapper->GetNumberOfColors() - 1;
    upperPixel = (unsigned char)(colors[upper_val].pixel);
    lowerPixel = (unsigned char)(colors[lower_val].pixel);
    }  

  inInc1 = -inInc1;

  // Loop through in regions pixels
  inPtr1 = inPtr;
  for (idx1 = inMin1; idx1 <= inMax1; idx1++)
    {
    inPtr0 = inPtr1;
    endPtr = inPtr0 + inInc0*(inMax0 - inMin0 + 1);
    if (visualClass == TrueColor)
      {
      while (inPtr0 != endPtr)
	{
	*ulOutPtr = 0;
	if (*inPtr0 <= lower)	  
	  {
	  *ulOutPtr = *ulOutPtr | ((rmask & (lowerPixel << 24)) >> rshift);
	  *ulOutPtr = *ulOutPtr | ((gmask & (lowerPixel << 24)) >> gshift);
	  *ulOutPtr = *ulOutPtr | ((bmask & (lowerPixel << 24)) >> bshift);
          ulOutPtr++;
	  }
	else if (*inPtr0 >= upper)
	  {
	  *ulOutPtr = *ulOutPtr | ((rmask & (upperPixel << 24)) >> rshift);
	  *ulOutPtr = *ulOutPtr | ((gmask & (upperPixel << 24)) >> gshift);
	  *ulOutPtr = *ulOutPtr | ((bmask & (upperPixel << 24)) >> bshift);
          ulOutPtr++;
	  }
	else
	  {
	  colorIdx = (int)((*inPtr0 + shift) * scale);
	  *ulOutPtr = *ulOutPtr | ((rmask & (colorIdx << 24)) >> rshift);
	  *ulOutPtr = *ulOutPtr | ((gmask & (colorIdx << 24)) >> gshift);
	  *ulOutPtr = *ulOutPtr | ((bmask & (colorIdx << 24)) >> bshift);
          ulOutPtr++;
	  }
	inPtr0 += inInc0;
	}
      }
    else if (visualClass == DirectColor)
      {
      while (inPtr0 != endPtr)
	{
        *ulOutPtr = 0;
	if (*inPtr0 <= lower) 
	  {
	  *ulOutPtr = *ulOutPtr | ((rmask & (lowerPixel << 24)) >> rshift);
	  *ulOutPtr = *ulOutPtr | ((gmask & (lowerPixel << 24)) >> gshift);
	  *ulOutPtr = *ulOutPtr | ((bmask & (lowerPixel << 24)) >> bshift);
          ulOutPtr++;
	  }
	else if (*inPtr0 >= upper)
	  {
	  *ulOutPtr = *ulOutPtr | ((rmask & (upperPixel << 24)) >> rshift);
	  *ulOutPtr = *ulOutPtr | ((gmask & (upperPixel << 24)) >> gshift);
	  *ulOutPtr = *ulOutPtr | ((bmask & (upperPixel << 24)) >> bshift);
          ulOutPtr++;
	  }
	else
	  {
	  colorIdx = (int)((*inPtr0 + shift) * scale);
	  *ulOutPtr = *ulOutPtr | ((rmask & (colorIdx << 24)) >> rshift);
	  *ulOutPtr = *ulOutPtr | ((gmask & (colorIdx << 24)) >> gshift);
	  *ulOutPtr = *ulOutPtr | ((bmask & (colorIdx << 24)) >> bshift);
          ulOutPtr++;
	  }
	inPtr0 += inInc0;
	}
      }
    else if (visualClass == PseudoColor)
      {
      while (inPtr0 != endPtr)
	{
	if (*inPtr0 <= lower) 
	  {
	  *outPtr++ = lowerPixel;
	  }
	else if (*inPtr0 >= upper)
	  {
	  *outPtr++ = upperPixel;
	  }
	else
	  {
	  colorIdx = (int)((*inPtr0 + shift) * scale);
	  *outPtr++ = (unsigned char)(colors[colorIdx].pixel);
	  }
	inPtr0 += inInc0;
	}
      }
    inPtr1 += inInc1;
    }

   delete[] colors;

}


//----------------------------------------------------------------------------
// Description:
// A templated function that handles color images. (only True Color 24 bit)
template <class T>
static void vtkXImageMapperRenderColor(vtkXImageMapper *mapper,
				       vtkViewport *viewport,
                                       vtkImageData *data, T *redPtr, int bpp,
                                       unsigned char *outPtr)
{
  int red, green, blue;
  T *redPtr0, *redPtr1;
  T *bluePtr0, *bluePtr1;
  T *greenPtr0, *greenPtr1;
  int inMin0, inMax0, inMin1, inMax1;
  int inInc0, inInc1;
  int idx0, idx1;
  float shift, scale;
  T *greenPtr;
  T *bluePtr;
  T lower, upper;
  unsigned char lower_val, upper_val;
  
  
  vtkWindow*  window = viewport->GetVTKWindow();
    
  shift = mapper->GetColorShift();
  scale = mapper->GetColorScale();

  int* tempExt = mapper->GetInput()->GetUpdateExtent();
  inMin0 = tempExt[0];
  inMax0 = tempExt[1];
  inMin1 = tempExt[2];
  inMax1 = tempExt[3];
  int* tempIncs = data->GetIncrements();
  inInc0 = tempIncs[0];
  inInc1 = tempIncs[1];

  if (bpp >= 2) greenPtr = redPtr + 1;
  else greenPtr = redPtr;
  if (bpp >= 3) bluePtr = redPtr + 2;
  else bluePtr = redPtr;

  unsigned long rmask = 0;
  unsigned long gmask = 0;
  unsigned long bmask = 0;

  mapper->GetXWindowColorMasks(window, &rmask, &gmask, &bmask);
  
  int rshift = 0;

  while ( ((rmask & 0x80000000) == 0) && (rshift < 32))
    {
    rmask = rmask << 1;
    rshift++;
    }

  int gshift = 0;

  while ( ((gmask & 0x80000000) == 0) && (gshift < 32))
    {
    gmask = gmask << 1;
    gshift++;
    }

  int bshift = 0;

  while ( ((bmask & 0x80000000) == 0) && (bshift < 32))
    {
    bmask = bmask << 1;
    bshift++;
    }

  unsigned long* ulOutPtr = (unsigned long*) outPtr;

  vtkXImageMapperClamps ( data, mapper->GetColorWindow(),
			  mapper->GetColorLevel(), 
			  lower, upper, lower_val, upper_val);

  inInc1 = -inInc1;

  // Loop through in regions pixels
  redPtr1 = redPtr;
  greenPtr1 = greenPtr;
  bluePtr1 = bluePtr;
  for (idx1 = inMin1; idx1 <= inMax1; idx1++)
    {
    redPtr0 = redPtr1;
    greenPtr0 = greenPtr1;
    bluePtr0 = bluePtr1;
    for (idx0 = inMin0; idx0 <= inMax0; idx0++)
      {
      if (*redPtr0 <= lower) red = lower_val;
      else if (*redPtr0 >= upper) red = upper_val;
      else red = (unsigned char)(((float)(*redPtr0) + shift) * scale);

      if (*greenPtr0 <= lower) green = lower_val;
      else if (*greenPtr0 >= upper) green = upper_val;
      else green = (unsigned char)(((float)(*greenPtr0) + shift) * scale);
  
      if (*bluePtr0 <= lower) blue = lower_val;
      else if (*bluePtr0 >= upper) blue = upper_val;
      else blue = (unsigned char)(((float)(*bluePtr0) + shift) * scale);

      *ulOutPtr = 0;
      *ulOutPtr = *ulOutPtr | ((rmask & (red << 24)) >> rshift);
      *ulOutPtr = *ulOutPtr | ((gmask & (green << 24)) >> gshift);
      *ulOutPtr = *ulOutPtr | ((bmask & (blue << 24)) >> bshift);
      ulOutPtr++;

      redPtr0 += inInc0;
      greenPtr0 += inInc0;
      bluePtr0 += inInc0;
      }
    redPtr1 += inInc1;
    greenPtr1 += inInc1;
    bluePtr1 += inInc1;
    }
}




void vtkXImageMapper::RenderData(vtkViewport* viewport, vtkImageData* data, vtkActor2D* actor)
{
  int width, height;
  int size;
  void *ptr0;

  vtkWindow*  window = viewport->GetVTKWindow();

  int visualDepth = this->GetXWindowDepth(window);
  Display* displayId = (Display*) window->GetGenericDisplayId();
  Visual visualId;  
  this->GetXWindowVisualId(window, &visualId);

  GC gc = (GC) window->GetGenericContext();
  if (gc == NULL) vtkErrorMacro(<<"Window returned NULL gc!");

  int* extent = this->GetInput()->GetUpdateExtent();
  width = (extent[1] - extent[0] + 1);
  height = (extent[3] - extent[2] + 1);

  size = width * height;

  // Allocate output data
  if (visualDepth >= 24)
  // if((visualClass == TrueColor) || (visualClass == DirectColor))
    // if (visualClass == TrueColor)
    {
    size *= 4;
    }

  // Only reallocate DataOut if the size is different than before
  if (size != this->DataOutSize)
    {
    if (this->DataOut)
      {
      delete[] this->DataOut;
      }
    this->DataOut = new unsigned char[size];
    if (!this->DataOut)
      {
      vtkErrorMacro(<< "Could not reallocate data out.");
      return;
      }
    this->DataOutSize = size;
    vtkDebugMacro(<<"vtkXImageMapper::RenderData - Reallocated DataOut");
    }

  int dim;
  dim = data->GetNumberOfScalarComponents();
  ptr0 = data->GetScalarPointer(extent[0], extent[3], extent[4]);

  if (dim > 1)  
    {
    // Call the appropriate templated function
    switch (data->GetScalarType())
      {
      case VTK_FLOAT:
        vtkDebugMacro(<<"vtkXImageMapper::RenderData - RenderColor, float");
        vtkXImageMapperRenderColor(this, viewport, data,
                           (float *)(ptr0), dim, this->DataOut);
        break;
      case VTK_INT:
        vtkDebugMacro(<<"vtkXImageMapper::RenderData - RenderColor, int");
        vtkXImageMapperRenderColor(this, viewport, data,
                           (int *)(ptr0), dim, this->DataOut);
        break;
      case VTK_SHORT:
        vtkDebugMacro(<<"vtkXImageMapper::RenderData - RenderColor, short");
        vtkXImageMapperRenderColor(this, viewport, data,
                           (short *)(ptr0),dim, this->DataOut);
        break;
      case VTK_UNSIGNED_SHORT:
        vtkDebugMacro(<<"vtkXImageMapper::RenderData - RenderColor, unsigned short");
        vtkXImageMapperRenderColor(this, viewport, data, (unsigned short *)(ptr0),
                            dim, this->DataOut);
        break;
      case VTK_UNSIGNED_CHAR:
        vtkDebugMacro(<<"vtkXImageMapper::RenderData - RenderColor, unsigned char ");
        vtkXImageMapperRenderColor(this, viewport, data, (unsigned char *)(ptr0),
                            dim, this->DataOut);
        break;
      }
    }
  else
    {
    // GrayScale images.
    //###   ptr0 = data->GetScalarPointer();
    // Call the appropriate templated function
    switch (data->GetScalarType())
      {
      case VTK_FLOAT:
        vtkXImageMapperRenderGray(this, viewport, data, 
				  (float *)(ptr0), this->DataOut);
        break;
      case VTK_INT:
        vtkXImageMapperRenderGray(this, viewport, data, 
				  (int *)(ptr0), this->DataOut);
        break;
      case VTK_SHORT:
        vtkXImageMapperRenderGray(this, viewport, data, 
				  (short *)(ptr0), this->DataOut);
        break;
      case VTK_UNSIGNED_SHORT:
        vtkXImageMapperRenderGray(this, viewport, data, 
				  (unsigned short *)(ptr0), this->DataOut);
        break;
      case VTK_UNSIGNED_CHAR:
        vtkXImageMapperRenderGray(this, viewport, data, 
				  (unsigned char *)(ptr0), this->DataOut);
        break;
      }  
    }

  vtkDebugMacro(<<"vtkXImageMapper::RenderData - Creating X image.");

  // Create the image pixmap.
  this->Image = XCreateImage(displayId, &visualId, visualDepth,
                             ZPixmap, 0, (char *)this->DataOut, width, height, 8,0);


  int* actPos = 
    actor->GetPositionCoordinate()->GetComputedLocalDisplayValue(viewport);
  
  // take into account adjustments
  actPos[0] += this->PositionAdjustment[0];
  actPos[1] -= this->PositionAdjustment[1];
  
  //float* actorScale = actor->GetScale();

  // The actor offset is the bottom left corner of the image
  // We need the upper left corner for displaying it on the screen,
  // so subtract the height of the image times it's Y scale.
  //actPos[1] = actPos[1] - (int)(actorScale[1]*height);
  actPos[1] = actPos[1] - height + 1;

  // Get the compositing mode for the actor
  int compositeMode = this->GetCompositingMode(actor);

  // Set the compositing mode in the graphics context (gc)
  XSetFunction(displayId, gc, compositeMode);

  // Put the image on the screen
  vtkDebugMacro(<<"vtkXImageMapper::RenderData - Putting X image on screen.");

  Drawable drawable = (Drawable) window->GetGenericDrawable();

  if (!drawable) vtkErrorMacro(<<"Window returned NULL drawable!");
  XPutImage(displayId, drawable, gc, this->Image, 0, 0,
            actPos[0], actPos[1], width, height);

  // Flush the X queue
  XFlush(displayId);
  XSync(displayId, False);
 
  // Deallocate the X image
  XFree(this->Image);

}




