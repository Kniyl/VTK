#!/usr/local/bin/python

from vtkpython import *
from WindowLevelInterface import *

# Shift and scale an image (in that order)
# This filter is usefull for converting to a lower precision data type.
# This script test the clamp overflow feature.



reader = vtkImageReader()
#reader.DebugOn()
reader.GetOutput().ReleaseDataFlagOff()
reader.SetDataByteOrderToLittleEndian()
reader.SetDataExtent(0,255,0,255,1,93)
reader.SetFilePrefix("../../../vtkdata/fullHead/headsq")
reader.SetDataMask(0x7fff)
#reader.GetOutput().SetMemoryLimit(30)

shiftScale = vtkImageShiftScale()
shiftScale.SetInput(reader.GetOutput())
shiftScale.SetShift(-1000.0)
shiftScale.SetScale(4.0)
shiftScale.SetOutputScalarTypeToUnsignedShort()
shiftScale.ClampOverflowOn()

shiftScale2 = vtkImageShiftScale()
shiftScale2.SetInput(shiftScale.GetOutput())
shiftScale2.SetShift(0)
shiftScale2.SetScale(2.0)

viewer = vtkImageViewer()
#viewer.DebugOn()
viewer.SetInput(shiftScale2.GetOutput())
#viewer.SetInput(reader.GetOutput())
viewer.SetColorWindow(1024)
viewer.SetColorLevel(512)

#make interface
WindowLevelInterface(viewer)

w = vtkPNMWriter()
w.SetFileName("D:/vtknew/vtk/graphics/examplesTcl/mace2.ppm")
w.SetInput(shiftScale.GetOutput())
#w.Write()
