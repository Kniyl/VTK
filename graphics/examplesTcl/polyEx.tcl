catch {load vtktcl}

# This example demonstrates the reading of point data and cell data
# simultaneously.

# get the interactor ui
source ../../examplesTcl/vtkInt.tcl

# create pipeline
#
vtkPolyDataReader reader
    reader SetFileName "../../../vtkdata/polyEx.vtk"
vtkPolyDataMapper mapper
    mapper SetInput [reader GetOutput]
    eval mapper SetScalarRange [[reader GetOutput] GetScalarRange]
vtkActor actor
    actor SetMapper mapper

# Create the RenderWindow, Renderer and both Actors
vtkRenderer ren1
vtkRenderWindow renWin
    renWin AddRenderer ren1
vtkRenderWindowInteractor iren
    iren SetRenderWindow renWin

ren1 AddActor actor
ren1 SetBackground 1 1 1
renWin SetSize 500 500

set cam1 [ren1 GetActiveCamera]
$cam1 SetClippingRange .348 17.43
$cam1 SetPosition 2.92 2.62 -0.836
$cam1 ComputeViewPlaneNormal
$cam1 SetViewUp -0.436 -0.067 -0.897

# render the image
#
iren SetUserMethod {wm deiconify .vtkInteract}
renWin Render

#renWin SetFileName "polyEx.tcl.ppm"
#renWin SaveImageAsPPM

# prevent the tk window from showing up then start the event loop
wm withdraw .


