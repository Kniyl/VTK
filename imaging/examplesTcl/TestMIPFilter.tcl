catch {load vtktcl}
# Simple viewer for images.


source vtkImageInclude.tcl

# Image pipeline

vtkImageReader reader
reader ReleaseDataFlagOff
reader SetDataByteOrderToLittleEndian
reader SetDataExtent 0 255 0 255 1 93
reader SetFilePrefix "../../../vtkdata/fullHead/headsq"
reader SetDataMask 0x7fff
reader DebugOn
#reader Update


vtkImageMIPFilter mip
mip SetInput [reader GetOutput]
mip MIPZOn
mip SetProjectionRange 0 92

vtkImageViewer viewer
viewer SetInput [mip GetOutput]
viewer SetColorWindow 3000
viewer SetColorLevel 1500
#viewer DebugOn
#viewer Render

#make interface
source WindowLevelInterface.tcl







