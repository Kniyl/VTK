#!/usr/local/bin/python

from vtkpython import *
from WindowLevelInterface import *

# Make an image larger by repeating the data.  Tile.

# Image pipeline
reader = vtkPNMReader()
reader.ReleaseDataFlagOff()
reader.SetFileName("../../../vtkdata/earth.ppm")

pad = vtkImageWrapPad()
pad.SetInput(reader.GetOutput())
pad.SetOutputWholeExtent(-100,100,0,250,0,0)

viewer = vtkImageViewer()
viewer.SetInput(pad.GetOutput())
viewer.SetZSlice(0)
viewer.SetColorWindow(255)
viewer.SetColorLevel(127)
viewer.GetActor2D().SetDisplayPosition(100,0)

#make interface
WindowLevelInterface(viewer)
