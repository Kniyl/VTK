#!/usr/local/bin/python

from libVTKCommonPython import *
from libVTKGraphicsPython import *
from libVTKImagingPython import *

#catch  load vtktcl 
# user interface command widget
#source ../../examplesTcl/vtkInt.tcl

# create a rendering window and renderer
ren = vtkRenderer()
renWin = vtkRenderWindow()
renWin.AddRenderer(ren)
renWin.SetSize(300,300)

# create an actor and give it cone geometry
cone = vtkConeSource()
cone.SetResolution(8)
coneMapper = vtkPolyDataMapper()
coneMapper.SetInput(cone.GetOutput())
coneActor = vtkActor()
coneActor.SetMapper(coneMapper)

# assign our actor to the renderer
ren.AddActor(coneActor)

# enable user interface interactor
renWin.Render()

w2if =  vtkWindowToImageFilter()
w2if.SetInput(renWin)

imgDiff = vtkImageDifference()
    
rtpnm = vtkPNMReader()
rtpnm.SetFileName("valid/Cone.tcl.ppm");

imgDiff.SetInput(w2if.GetOutput());
imgDiff.SetImage(rtpnm.GetOutput());
imgDiff.Update();

if imgDiff.GetThresholdedError() <= 10:
	print 'Python smoke test passed.'
else:
	print 'Python smoke test failed.'
