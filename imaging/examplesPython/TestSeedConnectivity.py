#!/usr/local/bin/python
import os
try:
  VTK_DATA = os.environ['VTK_DATA']
except KeyError:
  VTK_DATA = '../../../vtkdata/'


from vtkpython import *
from WindowLevelInterface import *

# A script to test the island removal filter.
# first the image is thresholded, then small islands are removed.


# Image pipeline

reader = vtkImageReader()
reader.ReleaseDataFlagOff()
reader.SetDataByteOrderToLittleEndian()
reader.SetDataExtent(0,255,0,255,1,93)
reader.SetFilePrefix(VTK_DATA + "/fullHead/headsq")
reader.SetDataMask(0x7fff)
#reader.DebugOn()

thresh = vtkImageThreshold()
thresh.SetInput(reader.GetOutput())
thresh.ThresholdByUpper(2000.0)
thresh.SetInValue(255)
thresh.SetOutValue(0)
thresh.ReleaseDataFlagOff()

cast = vtkImageCast()
cast.SetInput(thresh.GetOutput())
cast.SetOutputScalarType(VTK_UNSIGNED_CHAR)

connect = vtkImageSeedConnectivity()
connect.SetInput(cast.GetOutput())
connect.SetInputConnectValue(255)
connect.SetOutputConnectedValue(255)
connect.SetOutputUnconnectedValue(0)
connect.AddSeed(0,200,0)
#connect.DebugOn()

viewer = vtkImageViewer()
viewer.SetInput(connect.GetOutput())
viewer.SetZSlice(22)
viewer.SetColorWindow(255)
viewer.SetColorLevel(127.5)

# make interface
WindowLevelInterface(viewer)
