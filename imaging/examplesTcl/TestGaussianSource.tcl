catch {load vtktcl}
if { [catch {set VTK_TCL $env(VTK_TCL)}] != 0} { set VTK_TCL "../../examplesTcl" }
if { [catch {set VTK_DATA $env(VTK_DATA)}] != 0} { set VTK_DATA "../../../vtkdata" }

# A script to test the GaussianSource

source vtkImageInclude.tcl

# Image pipeline

vtkImageGaussianSource gauss
gauss SetWholeExtent 0 225 0 225 0 20
gauss SetCenter 100 100 10
gauss SetStandardDeviation 100.0
gauss SetMaximum 255.0
gauss ReleaseDataFlagOff

vtkImageViewer viewer
viewer SetInput [gauss GetOutput]
viewer SetZSlice 10
viewer SetColorWindow 255
viewer SetColorLevel 127.5
#viewer DebugOn

# make interface
source WindowLevelInterface.tcl







