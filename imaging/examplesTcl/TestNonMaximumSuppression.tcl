# This script is for testing the 3d NonMaximumSuppressionFilter.
# The filter works exclusively on the output of the gradient filter.
# The effect is to pick the peaks of the gradient creating thin surfaces.

catch {load vtktcl}
source vtkImageInclude.tcl



# Image pipeline
vtkImageReader reader
reader SetDataByteOrderToLittleEndian
reader SetDataExtent 0 255 0 255 1 93
reader SetFilePrefix "../../../vtkdata/fullHead/headsq"
reader SetDataMask 0x7fff
#reader DebugOn

vtkImageGradient gradient
gradient SetFilteredAxes $VTK_IMAGE_X_AXIS $VTK_IMAGE_Y_AXIS $VTK_IMAGE_Z_AXIS
gradient SetInput [reader GetOutput]
gradient ReleaseDataFlagOff

vtkImageMagnitude magnitude
magnitude SetInput [gradient GetOutput]

vtkImageNonMaximumSuppression suppress
suppress SetVectorInput [gradient GetOutput]
suppress SetMagnitudeInput [magnitude GetOutput]
suppress SetNumberOfFilteredAxes 3
suppress ReleaseDataFlagOff

vtkImageViewer viewer
viewer SetInput [suppress GetOutput]
viewer SetCoordinate2 22
viewer SetColorWindow 1000
viewer SetColorLevel 500


# make interface
source WindowLevelInterface.tcl






