catch {load vtktcl}
if { [catch {set VTK_TCL $env(VTK_TCL)}] != 0} { set VTK_TCL "../../examplesTcl" }
if { [catch {set VTK_DATA $env(VTK_DATA)}] != 0} { set VTK_DATA "../../../vtkdata" }

# This script calculates the luminanace of an image

vtkImageWindow imgWin


# Image pipeline

vtkBMPReader image
  image SetFileName "$VTK_DATA/beach.bmp"

vtkImageShrink3D shrink
shrink SetInput [image GetOutput]
shrink SetShrinkFactors 4 4 1

set operators "\
ByPass \
And \
Nand \
Xor \ 
Or \ 
Nor"

foreach operator $operators {
    vtkImageMaskBits operator${operator}
      operator${operator} SetInput [shrink GetOutput]
      if { $operator != "ByPass" } {
        operator${operator} SetOperationTo${operator}
      } else {
        operator${operator} BypassOn
      }     
      operator${operator} SetMasks 255 255 0
    vtkImageMapper mapper${operator}
      mapper${operator} SetInput [operator${operator} GetOutput]
      mapper${operator} SetColorWindow 255
      mapper${operator} SetColorLevel 127.5
    vtkActor2D actor${operator}
      actor${operator} SetMapper mapper${operator}
    vtkImager imager${operator}
      imager${operator} AddActor2D actor${operator}
    imgWin AddImager imager${operator}
}

set column 1
set row 1
set deltaX [expr 1.0/3.0]
set deltaY [expr 1.0/2.0]

foreach operator $operators {
    imager${operator} SetViewport [expr ($column - 1) * $deltaX] [expr ($row - 1) * $deltaY] [expr $column * $deltaX] [expr $row * $deltaY]
    incr column
    if { $column > 3 } {set column 1; incr row}
}

imgWin SetSize 384 256
imgWin Render
imgWin SetFileName TestAllMaskBits.tcl.ppm
#imgWin SaveImageAsPPM

wm withdraw .






