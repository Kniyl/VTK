package require vtktcl

# This script calculates the luminanace of an image

vtkImageWindow imgWin


# Image pipeline

vtkTIFFReader image1
  image1 SetFileName "$VTK_DATA_ROOT/Data/beach.tif"

vtkImageShrink3D shrink
shrink SetInput [image1 GetOutput]
shrink SetShrinkFactors 2 2 1

set operators "\
ByPass \
And \
Nand \
Xor \ 
Or \ 
Nor"

foreach operator $operators {
    if { $operator != "ByPass" } {    
	vtkImageMaskBits operator${operator}
	operator${operator} SetInput [shrink GetOutput]
        operator${operator} SetOperationTo${operator}
	operator${operator} SetMasks 255 255 0
    } 
    
    vtkImageMapper mapper${operator}
    if { $operator != "ByPass" } {    	
	mapper${operator} SetInput [operator${operator} GetOutput]
    } else {
	mapper${operator} SetInput [shrink GetOutput]
    }
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

wm withdraw .






