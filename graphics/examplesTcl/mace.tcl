catch {load vtktcl}
#
# this is a tcl version of the Mace example
# get the interactor ui
source ../../examplesTcl/vtkInt.tcl

# Create the RenderWindow, Renderer and both Actors
#
vtkRenderer ren1
vtkRenderWindow renWin
    renWin AddRenderer ren1
    renWin SetWindowName "vtk - Mace"
vtkRenderWindowInteractor iren
    iren SetRenderWindow renWin

# create a sphere source and actor
#
vtkSphereSource sphere

vtkPolyDataMapper   sphereMapper
    sphereMapper SetInput [sphere GetOutput]
    sphereMapper GlobalImmediateModeRenderingOn
vtkLODActor sphereActor
    sphereActor SetMapper sphereMapper

# create the spikes using a cone source and the sphere source
#
vtkConeSource cone
vtkGlyph3D glyph
    glyph SetInput [sphere GetOutput]
    glyph SetSource [cone GetOutput]
    glyph SetVectorModeToUseNormal
    glyph SetScaleModeToScaleByVector
    glyph SetScaleFactor 0.25
vtkPolyDataMapper spikeMapper
    spikeMapper SetInput [glyph GetOutput]
vtkLODActor spikeActor
    spikeActor SetMapper spikeMapper

# Add the actors to the renderer, set the background and size
#
ren1 AddActor sphereActor
ren1 AddActor spikeActor
ren1 SetBackground 0.1 0.2 0.4
renWin SetSize 300 300

# render the image
#
iren SetUserMethod {wm deiconify .vtkInteract}

set cam1 [ren1 GetActiveCamera]
$cam1 Zoom 1.4
iren Initialize
#renWin SetFileName "mace.tcl.ppm"
#renWin SaveImageAsPPM

proc TkCheckAbort {} {
  set foo [renWin GetEventPending]
  if {$foo != 0} {renWin SetAbortRender 1}
}
renWin SetAbortCheckMethod {TkCheckAbort}

# prevent the tk window from showing up then start the event loop
wm withdraw .


