catch {load vtktcl}
if { [catch {set VTK_TCL $env(VTK_TCL)}] != 0} { set VTK_TCL "../../examplesTcl" }
if { [catch {set VTK_DATA $env(VTK_DATA)}] != 0} { set VTK_DATA "../../../vtkdata" }

# Derived from Cursor3D.  This script increases the coverage of the
# vtkImageInplaceFilter super class.


source vtkImageInclude.tcl

# global values
set CURSOR_X 20
set CURSOR_Y 20
set CURSOR_Z 20

set IMAGE_MAG_X 2
set IMAGE_MAG_Y 2
set IMAGE_MAG_Z 1



# pipeline stuff
vtkSLCReader reader
    reader SetFileName "$VTK_DATA/poship.slc"

# make the image a little biger
vtkImageMagnify magnify1
  magnify1 SetInput [reader GetOutput]
  magnify1 SetMagnificationFactors $IMAGE_MAG_X $IMAGE_MAG_Y $IMAGE_MAG_Z
  magnify1 ReleaseDataFlagOn

vtkImageMagnify magnify2
  magnify2 SetInput [reader GetOutput]
  magnify2 SetMagnificationFactors $IMAGE_MAG_X $IMAGE_MAG_Y $IMAGE_MAG_Z
  magnify2 ReleaseDataFlagOn

# a filter that does in place processing (magnify ReleaseDataFlagOn)
vtkImageCursor3D cursor
cursor SetInput [magnify1 GetOutput]
cursor SetCursorPosition [expr $CURSOR_X * $IMAGE_MAG_X] \
    [expr $CURSOR_Y * $IMAGE_MAG_Y] [expr $CURSOR_Z * $IMAGE_MAG_Z]
cursor SetCursorValue 255
cursor SetCursorRadius [expr 50 * $IMAGE_MAG_X]
# stream to increase coverage of in place filter.

# put thge two together in one image
vtkImageAppend imageAppend
imageAppend SetAppendAxis 0
imageAppend AddInput [magnify2 GetOutput]
imageAppend AddInput [cursor GetOutput]

vtkImageViewer viewer
viewer SetInput [imageAppend GetOutput]
viewer SetZSlice [expr $CURSOR_Z * $IMAGE_MAG_Z]
viewer SetColorWindow 200
viewer SetColorLevel 80
#viewer DebugOn
viewer Render

viewer SetPosition 50 50

#make interface
source WindowLevelInterface.tcl

