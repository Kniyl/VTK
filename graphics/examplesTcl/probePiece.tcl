catch {load vtktcl}
if { [catch {set VTK_TCL $env(VTK_TCL)}] != 0} { set VTK_TCL "../../examplesTcl" }
if { [catch {set VTK_DATA $env(VTK_DATA)}] != 0} { set VTK_DATA "../../../vtkdata" }

# get the interactor ui
source $VTK_TCL/vtkInt.tcl
source $VTK_TCL/colors.tcl

# Create the RenderWindow, Renderer and both Actors
#
vtkRenderer ren1
vtkRenderWindow renWin
    renWin AddRenderer ren1
vtkRenderWindowInteractor iren
    iren SetRenderWindow renWin

# create pipeline
#
vtkImageReader reader
reader SetDataByteOrderToLittleEndian
reader SetDataExtent 0 127 0 127 1 93
reader SetFilePrefix "$VTK_DATA/headsq/half"
reader SetDataMask 0x7fff
reader SetDataSpacing 1.6 1.6 1.5

vtkContourFilter iso
iso SetInput [reader GetOutput]
iso SetValue 0 1150

vtkPolyDataMapper isoMapper
isoMapper SetInput [iso GetOutput]
isoMapper ScalarVisibilityOff
isoMapper SetPiece 2
isoMapper SetNumberOfPieces 4

vtkActor isoActor
isoActor SetMapper isoMapper
eval [isoActor GetProperty] SetColor $antique_white


# Add the actors to the renderer, set the background and size
#
ren1 AddActor isoActor
ren1 SetBackground 0.2 0.3 0.4
renWin SetSize 450 450
[ren1 GetActiveCamera] Elevation 90
[ren1 GetActiveCamera] SetViewUp 0 0 -1
ren1 ResetCameraClippingRange
renWin Render

# render the image
#
iren SetUserMethod {wm deiconify .vtkInteract}
#renWin SetFileName "valid/mcubes.tcl.ppm"
#renWin SaveImageAsPPM

# prevent the tk window from showing up then start the event loop
wm withdraw .


