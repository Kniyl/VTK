#!/usr/local/bin/python

from libVTKCommonPython import *
from libVTKGraphicsPython import *

#catch  load vtktcl 
# Created oriented text
#source ../../examplesTcl/vtkInt.tcl

# pipeline
axes = vtkAxes()
axes.SetOrigin(0,0,0)
axesMapper = vtkPolyDataMapper()
axesMapper.SetInput(axes.GetOutput())
axesActor = vtkActor()
axesActor.SetMapper(axesMapper)

atext = vtkVectorText()
atext.SetText("Origin")
textMapper = vtkPolyDataMapper()
textMapper.SetInput(atext.GetOutput())
textActor = vtkFollower()
textActor.SetMapper(textMapper)
textActor.SetScale(0.2,0.2,0.2)
textActor.AddPosition(0,-0.1,0)

# create graphics stuff
ren = vtkRenderer()
renWin = vtkRenderWindow()
renWin.AddRenderer(ren)
renWin.SetSize(400,400)
iren = vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin)

ren.AddActor(axesActor)
ren.AddActor(textActor)
ren.GetActiveCamera().Zoom(1.6)
renWin.Render()
textActor.SetCamera(ren.GetActiveCamera())

iren.Initialize()

renWin.SetFileName("textOrigin.tcl.ppm")
#renWin SaveImageAsPPM

#wm withdraw .
iren.Start()
