# Tst the OpenClose3D filter.

catch {load vtktcl}
source vtkImageInclude.tcl

# Image pipeline

vtkImageReader reader
reader SetDataByteOrderToLittleEndian
reader SetDataExtent 0 255 0 255 1 93
reader SetFilePrefix "../../../vtkdata/fullHead/headsq"
reader SetDataMask 0x7fff
#reader DebugOn

vtkImageThreshold thresh
thresh SetInput [reader GetOutput]
thresh SetOutputScalarTypeToUnsignedChar
thresh ThresholdByUpper 2000.0
thresh SetInValue 255
thresh SetOutValue 0
thresh ReleaseDataFlagOff

vtkImageOpenClose3D my_close
my_close SetInput [thresh GetOutput]
my_close SetOpenValue 0
my_close SetCloseValue 255
my_close SetKernelSize 5 5 3
my_close ReleaseDataFlagOff
# for coverage (we could compare results to see if they are correct).
my_close DebugOn
my_close DebugOff
my_close GetOutput
my_close GetCloseValue
my_close GetOpenValue
#my_close SetProgressMethod {set pro [my_close GetProgress]; puts "Completed $pro"; flush stdout}

vtkImageViewer viewer
viewer SetInput [my_close GetOutput]
viewer SetZSlice 2
viewer SetColorWindow 255
viewer SetColorLevel 127.5


# make interface
source WindowLevelInterface.tcl







