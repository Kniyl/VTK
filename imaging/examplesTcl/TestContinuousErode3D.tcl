# Min of sphere around every pixel.

catch {load vtktcl}
source vtkImageInclude.tcl

# Image pipeline

vtkImageReader reader
reader ReleaseDataFlagOff
reader SetDataByteOrderToLittleEndian
reader SetDataExtent 0 255 0 255 1 93
reader SetFilePrefix "../../../vtkdata/fullHead/headsq"
reader SetDataMask 0x7fff
#reader DebugOn

vtkImageContinuousErode3D erode
erode SetInput [reader GetOutput]
erode SetKernelSize 5 5 5

vtkImageViewer viewer
viewer SetInput [erode GetOutput]
viewer SetZSlice 22
viewer SetColorWindow 2000
viewer SetColorLevel 1000
#viewer DebugOn


source WindowLevelInterface.tcl


