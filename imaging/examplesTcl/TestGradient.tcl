catch {load vtktcl}
if { [catch {set VTK_TCL $env(VTK_TCL)}] != 0} { set VTK_TCL "../../examplesTcl" }
if { [catch {set VTK_DATA $env(VTK_DATA)}] != 0} { set VTK_DATA "../../../vtkdata" }

# This script is for testing the 2d Gradient filter.
# It only displays the first component (0) which contains
# the magnitude of the gradient.



source vtkImageInclude.tcl


# Image pipeline

vtkImageReader reader
#reader DebugOn
reader SetDataByteOrderToLittleEndian
reader SetDataExtent 0 255 0 255 1 93
reader SetFilePrefix "$VTK_DATA/fullHead/headsq"
reader SetDataMask 0x7fff

vtkImageGradient gradient
gradient SetInput [reader GetOutput]
gradient SetDimensionality 3

vtkImageViewer viewer
#viewer DebugOn
viewer SetInput [gradient GetOutput]
viewer SetZSlice 22
viewer SetColorWindow 400
viewer SetColorLevel 0


#make interface
source WindowLevelInterface.tcl







