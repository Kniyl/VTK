#!/usr/local/bin/python

from libVTKCommonPython import *
from libVTKGraphicsPython import *

#catch  load vtktcl 
# get the interactor ui
#source ../../examplesTcl/vtkInt.tcl

# create tensor ellipsoids
# Create the RenderWindow, Renderer and interactive renderer
#
ren = vtkRenderer()
renWin = vtkRenderWindow()
renWin.AddRenderer(ren)
iren = vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin)

#
# Create tensor ellipsoids
#
# generate tensors
ptLoad = vtkPointLoad()
ptLoad.SetLoadValue(100.0)
ptLoad.SetSampleDimensions(6,6,6)
ptLoad.ComputeEffectiveStressOn()
ptLoad.SetModelBounds(-10,10,-10,10,-10,10)

# extract plane of data
plane = vtkStructuredPointsGeometryFilter()
plane.SetInput(ptLoad.GetOutput())
plane.SetExtent(2,2,0,99,0,99)

# Generate ellipsoids
axes = vtkAxes()
axes.SetScaleFactor(0.5)
tubeAxes = vtkTubeFilter()
tubeAxes.SetInput(axes.GetOutput())
tubeAxes.SetRadius(0.1)
tubeAxes.SetNumberOfSides(6)
ellipsoids = vtkTensorGlyph()
ellipsoids.SetInput(ptLoad.GetOutput())
ellipsoids.SetSource(axes.GetOutput())
ellipsoids.SetScaleFactor(10)
ellipsoids.ClampScalingOn()
  
# Map contour
lut = vtkLogLookupTable()
lut.SetHueRange(.6667,0.0)
ellipMapper = vtkPolyDataMapper()
ellipMapper.SetInput(ellipsoids.GetOutput())
ellipMapper.SetLookupTable(lut)
plane.Update() #force.update.for.scalar.range()
ellipMapper.SetScalarRange(plane.GetOutput().GetScalarRange())

ellipActor = vtkActor()
ellipActor.SetMapper(ellipMapper)
ellipActor.GetProperty().SetAmbient(1)
ellipActor.GetProperty().SetDiffuse(0)
#
# Create outline around data
#
outline = vtkOutlineFilter()
outline.SetInput(ptLoad.GetOutput())

outlineMapper = vtkPolyDataMapper()
outlineMapper.SetInput(outline.GetOutput())

outlineActor = vtkActor()
outlineActor.SetMapper(outlineMapper)
outlineActor.GetProperty().SetColor(0,0,0)

#
# Create cone indicating application of load
#
coneSrc = vtkConeSource()
coneSrc.SetRadius(.5)
coneSrc.SetHeight(2)
coneMap = vtkPolyDataMapper()
coneMap.SetInput(coneSrc.GetOutput())
coneActor = vtkActor()
coneActor.SetMapper(coneMap)
coneActor.SetPosition(0,0,11)
coneActor.RotateY(90)
coneActor.GetProperty().SetColor(1,0,0)

camera = vtkCamera()
camera.SetFocalPoint(0.113766,-1.13665,-1.01919)
camera.SetPosition(-29.4886,-63.1488,26.5807)
camera.ComputeViewPlaneNormal()
camera.SetViewAngle(24.4617)
camera.SetViewUp(0.17138,0.331163,0.927879)

ren.AddActor(ellipActor)
ren.AddActor(outlineActor)
ren.AddActor(coneActor)
ren.SetBackground(1.0,1.0,1.0)
ren.SetActiveCamera(camera)

renWin.SetSize(500,500)
renWin.Render()

#renWin SetFileName TenAxes.tcl.ppm
#renWin SaveImageAsPPM


# prevent the tk window from showing up then start the event loop
#wm withdraw .




iren.Start()
