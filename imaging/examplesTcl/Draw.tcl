# Use the CanvasSource to draw in gray scale.
catch {load vtktcl}
if { [catch {set VTK_TCL $env(VTK_TCL)}] != 0} { set VTK_TCL "../../examplesTcl" }
if { [catch {set VTK_DATA $env(VTK_DATA)}] != 0} { set VTK_DATA "$VTK_DATA" }


source vtkImageInclude.tcl

vtkImageCanvasSource2D canvas
canvas SetScalarType $VTK_UNSIGNED_CHAR
canvas SetExtent 0 511 0 511 0 0
canvas SetDrawColor 66
canvas FillBox 0 511 0 511
canvas SetDrawColor 132
canvas FillBox 32 511 100 500
canvas SetDrawColor 33
canvas FillTube 500 20 30 400 5
canvas SetDrawColor 255
canvas DrawSegment 10 20 90 510
canvas SetDrawColor 100
canvas DrawSegment 510 90 10 20

# Check segment clipping
canvas SetDrawColor 66
canvas DrawSegment -10 30 30 -10
canvas DrawSegment -10 481 30 521
canvas DrawSegment 481 -10 521 30
canvas DrawSegment 481 521 521 481

# Check Filling a triangle
canvas SetDrawColor 70
canvas FillTriangle 100 100  300 150  150 300

# Check drawing a circle
canvas SetDrawColor 84
canvas DrawCircle 350 350  200.0

# Check drawing a point
canvas SetDrawColor 250
canvas DrawPoint 350 350

# Test filling functionality
canvas SetDrawColor 48
canvas DrawCircle 450 350 80.0
canvas SetDrawColor 255
canvas FillPixel 450 350


vtkImageViewer viewer
viewer SetInput [canvas GetOutput]
viewer SetColorWindow 256
viewer SetColorLevel 127.5


# make interface
source WindowLevelInterface.tcl






