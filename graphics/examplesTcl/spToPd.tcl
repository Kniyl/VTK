#
# structured points to geometry
#

catch {load vtktcl}
if { [catch {set VTK_TCL $env(VTK_TCL)}] != 0} { set VTK_TCL "../../examplesTcl" }
if { [catch {set VTK_DATA $env(VTK_DATA)}] != 0} { set VTK_DATA "../../../vtkdata" }

source $VTK_TCL/vtkInt.tcl

# Create the RenderWindow, Renderer and both Actors
#
vtkRenderer ren1
vtkRenderWindow renWin
    renWin AddRenderer ren1
vtkRenderWindowInteractor iren
    iren SetRenderWindow renWin

# read in some structured points
#
vtkPNMReader reader
  reader SetFileName $VTK_DATA/B.pgm

vtkStructuredPointsGeometryFilter geometry
  geometry SetInput [reader GetOutput]
  geometry SetExtent 0 10000 0 10000 0 0

vtkWarpScalar warp
  warp SetInput [geometry GetOutput]
  warp SetScaleFactor -.1

vtkDataSetMapper mapper
  mapper SetInput [warp GetOutput]
  mapper SetScalarRange 0 255
  mapper ImmediateModeRenderingOff

vtkActor actor
  actor SetMapper mapper


# Add the actors to the renderer, set the background and size
#
ren1 AddActor actor
[ren1 GetActiveCamera] Azimuth 20
[ren1 GetActiveCamera] Elevation 30
ren1 SetBackground 0.1 0.2 0.4
renWin SetSize 640 480

# render the image
#
iren SetUserMethod {wm deiconify .vtkInteract}
set cam1 [ren1 GetActiveCamera]
$cam1 Zoom 1.4
ren1 ResetCameraClippingRange
renWin Render
#renWin SetFileName "spToPd.tcl.ppm"
#renWin SaveImageAsPPM

# prevent the tk window from showing up then start the event loop
wm withdraw .


