catch {load vtktcl}
if { [catch {set VTK_TCL $env(VTK_TCL)}] != 0} { set VTK_TCL "../../examplesTcl" }
if { [catch {set VTK_DATA $env(VTK_DATA)}] != 0} { set VTK_DATA "../../../vtkdata" }

# This script shows the result of an ideal highpass filter in frequency space.

source vtkImageInclude.tcl

# Image pipeline

vtkImageReader reader
reader SetDataByteOrderToLittleEndian
reader SetDataExtent 0 255 0 255 1 93
reader SetFilePrefix "$VTK_DATA/fullHead/headsq"
reader SetDataMask 0x7fff
#reader DebugOn

vtkImageFFT fft
fft SetDimensionality 2
fft SetInput [reader GetOutput] 
#fft DebugOn

vtkImageButterworthHighPass highPass
highPass SetInput [fft GetOutput]
highPass SetOrder 2
highPass SetXCutOff 0.2
highPass SetYCutOff 0.1
highPass ReleaseDataFlagOff
#highPass DebugOn

vtkImageViewer viewer
viewer SetInput [highPass GetOutput]
viewer SetZSlice 22
viewer SetColorWindow 10000
viewer SetColorLevel 5000
#viewer DebugOn


# make interface
source WindowLevelInterface.tcl


