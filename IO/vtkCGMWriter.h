/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkCGMWriter.h
  Language:  C++
  Date:      $Date$
  Version:   $Revision$
  Credit:    The origin of much of this code was from the cd package
             written by G. Edward Johnson at the National Institute 
             of Standards and Technology (US).


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
// .NAME vtkCGMWriter - write polygonal data as a CGM file
// .SECTION Description
// vtkCGMWriter writes CGM (Computer Graphics Metafile) output. CGM is a 2D
// graphics vector format typically used by large plotters. This writer can 
// handle vertices, lines, polygons, and triangle strips in any combination. 
// Colors are specified either 1) from cell scalars (assumed to be RGB or 
// RGBA color specification), 2) from a specified color; or 3) randomly 
// assigned colors.
//
// Note: During output of the polygonal data, stringle strips are converted
// to triangles, and polylines to lines. Also, due to limitations in the CGM 
// color model, only 256 colors are available to the color palette.

// .SECTION Caveats
// The class vtkImageToPolyDataFilter is convenient for converting a raster
// image into polygons (and color map) suitable for plotting with CGM.

// .SECTION See Also
// vtkPolyDataWriter vtkPointDataToCellData


#ifndef __vtkCGMWriter_h
#define __vtkCGMWriter_h

#include "vtkPolyDataWriter.h"
#include "vtkViewport.h"

#define VTK_COLOR_MODE_DEFAULT 0
#define VTK_COLOR_MODE_SPECIFIED_COLOR 1
#define VTK_COLOR_MODE_RANDOM_COLORS 2

class VTK_IO_EXPORT vtkCGMWriter : public vtkPolyDataWriter
{
public:
  // Description:
  // Instantiate with no viewport defined and sorting on. The default
  // resolution is 10,000, and the color mode is set to default.
  static vtkCGMWriter *New() {return new vtkCGMWriter;};

  vtkTypeMacro(vtkCGMWriter,vtkPolyDataWriter);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Specify a vtkViewport object to be used to transform the vtkPolyData
  // points into 2D coordinates. By default (no vtkViewport specified), the 
  // point coordinates are generated by ignoring the z values. If a viewport
  // is defined, then the points are transformed into viewport coordinates.
  vtkSetObjectMacro(Viewport, vtkViewport);
  vtkGetObjectMacro(Viewport, vtkViewport);

  // Description:
  // Turn on/off the sorting of the cells via depth. If enabled, polygonal
  // cells will be sorted from back to front, i.e., a Painter's algorithm
  // sort.
  vtkSetMacro(Sort,int);
  vtkGetMacro(Sort,int);

  // Description:
  // Specify the resolution of the CGM file. This number is used to integerize
  // the maximum coordinate range of the plot file.
  vtkSetClampMacro(Resolution, int, 100, VTK_LARGE_INTEGER);
  vtkGetMacro(Resolution, int);

  // Description:
  // Control how output polydata is colored. By default (ColorModeToDefault),
  // if per cell colors are defined (unsigned chars of 1-4 components), then
  // the cells are colored with these values. (If point colors are defined
  // and cell colors are not, you can use vtkPointDataToCellData to convert
  // the point colors to cell colors.) Otherwise, by default, the cells are
  // set to the specified color. If ColorModeToSpecifiedColor is set, then
  // the primitives will all be set to this color. If ColorModeToRandomColors
  // is set, each cell will be randomly assigned a color.
  vtkSetMacro(ColorMode,int);
  vtkGetMacro(ColorMode,int);
  void SetColorModeToDefault() {
    this->SetColorMode(VTK_COLOR_MODE_DEFAULT);};
  void SetColorModeToSpecifiedColor() {
    this->SetColorMode(VTK_COLOR_MODE_SPECIFIED_COLOR);};
  void SetColorModeToRandomColors() {
    this->SetColorMode(VTK_COLOR_MODE_RANDOM_COLORS);};

  // Description:
  // Set/Get the specified color to color the polydata cells. This
  // color is only used when the color mode is set to 
  // ColorModeToSpecifiedColor, or ColorModeToDefault is set and no
  // cell colors are specified. The specified color is specified as RGB 
  // values ranging from (0,1). (Note: CGM will map this color to the
  // closest color it supports.)
  vtkSetVector3Macro(SpecifiedColor,float);
  vtkGetVectorMacro(SpecifiedColor,float,3);

protected:
  vtkCGMWriter();
  ~vtkCGMWriter();
  vtkCGMWriter(const vtkCGMWriter&);
  void operator=(const vtkCGMWriter&);
  void WriteData();

  vtkViewport *Viewport;
  int         ColorMode;
  float       SpecifiedColor[3];
  int         Resolution;
  int         Sort;
  
};

#endif

