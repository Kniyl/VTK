catch {load vtktcl}
# This script shows the result of an ideal lowpass filter in frequency space.

source vtkImageInclude.tcl

# Image pipeline

vtkImageReader reader
reader SetDataByteOrderToLittleEndian
reader SetDataExtent 0 255 0 255 1 93
reader SetFilePrefix "../../../vtkdata/fullHead/headsq"
reader SetDataMask 0x7fff
#reader DebugOn

vtkImageFFT fft
fft SetFilteredAxes $VTK_IMAGE_X_AXIS $VTK_IMAGE_Y_AXIS
fft SetInput [reader GetOutput]
#fft DebugOn

vtkImageButterworthLowPass lowPass
lowPass SetInput [fft GetOutput]
lowPass SetOrder 2
lowPass SetXCutOff 0.2
lowPass SetYCutOff 0.1
lowPass ReleaseDataFlagOff
#lowPass DebugOn

vtkImageViewer viewer
viewer SetInput [lowPass GetOutput]
viewer SetZSlice 22
viewer SetColorWindow 10000
viewer SetColorLevel 5000
#viewer DebugOn


# make interface
source WindowLevelInterface.tcl


