catch {load vtktcl}
# A script to test DilationErode filter
# First the image is thresholded.
# It is the dilated with a spher of radius 5.

source vtkImageInclude.tcl

# Image pipeline

vtkImageReader reader
reader SetDataByteOrderToLittleEndian
reader SetDataExtent 0 255 0 255 1 93
reader SetFilePrefix "../../../vtkdata/fullHead/headsq"
reader SetDataMask 0x7fff
reader SetOutputScalarType $VTK_SHORT
reader DebugOn

vtkImageThreshold thresh
thresh SetInput [reader GetOutput]
thresh SetOutputScalarType $VTK_UNSIGNED_CHAR
thresh ThresholdByUpper 2000.0
thresh SetInValue 255
thresh SetOutValue 0

vtkImageDilateErode3D dilate
dilate SetInput [thresh GetOutput]
dilate SetDilateValue 255
dilate SetErodeValue 0
dilate SetKernelSize 5 5 5
dilate ReleaseDataFlagOff

vtkImageViewer viewer
viewer SetInput [dilate GetOutput]
viewer SetZSlice 22
viewer SetColorWindow 255
viewer SetColorLevel 128
#viewer DebugOn


# make interface
source WindowLevelInterface.tcl
