# This script is for testing the 3d NonMaximumSuppressionFilter.
# The filter works exclusively on the output of the gradient filter.
# The effect is to pick the peaks of the gradient creating thin surfaces.

catch {load vtktcl}
if { [catch {set VTK_TCL $env(VTK_TCL)}] != 0} { set VTK_TCL "../../examplesTcl" }
if { [catch {set VTK_DATA $env(VTK_DATA)}] != 0} { set VTK_DATA "../../../vtkdata" }

source vtkImageInclude.tcl



# Image pipeline
vtkImageReader reader
reader SetDataByteOrderToLittleEndian
reader SetDataExtent 0 255 0 255 1 93
reader SetFilePrefix "$VTK_DATA/fullHead/headsq"
reader SetDataMask 0x7fff
#reader DebugOn

vtkImageGradient gradient
gradient SetDimensionality 2
gradient SetInput [reader GetOutput]

vtkImageMagnitude magnitude
magnitude SetInput [gradient GetOutput]

vtkImageNonMaximumSuppression suppress
suppress SetVectorInput [gradient GetOutput]
suppress SetMagnitudeInput [magnitude GetOutput]
suppress SetDimensionality 2

vtkImageViewer viewer
viewer SetInput [suppress GetOutput]
viewer SetZSlice 22
viewer SetColorWindow 1000
viewer SetColorLevel 500


# make interface
source WindowLevelInterface.tcl






