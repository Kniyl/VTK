# This Script test the euclidean to polar by coverting 2D vectors 
# from a gradient into polar, which is converted into HSV, and then to RGB.
catch {load vtktcl}
if { [catch {set VTK_TCL $env(VTK_TCL)}] != 0} { set VTK_TCL "../../examplesTcl" }
if { [catch {set VTK_DATA $env(VTK_DATA)}] != 0} { set VTK_DATA "../../../vtkdata" }


source vtkImageInclude.tcl


# Image pipeline

vtkImageGaussianSource gauss
gauss SetWholeExtent 0 255 0 255 0 44
gauss SetCenter 127 127 22 
gauss SetStandardDeviation 50.0
gauss SetMaximum 10000.0

vtkImageGradient gradient
gradient SetInput [gauss GetOutput]
gradient SetDimensionality 2

vtkImageEuclideanToPolar polar
polar SetInput [gradient GetOutput]
polar SetThetaMaximum 255

vtkImageViewer viewer
#viewer DebugOn
viewer SetInput [polar GetOutput]
viewer SetZSlice 22
viewer SetColorWindow 255
viewer SetColorLevel 127.5


#make interface
source WindowLevelInterface.tcl







