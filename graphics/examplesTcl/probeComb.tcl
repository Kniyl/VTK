catch {load vtktcl}
if { [catch {set VTK_TCL $env(VTK_TCL)}] != 0} { set VTK_TCL "../../examplesTcl" }
if { [catch {set VTK_DATA $env(VTK_DATA)}] != 0} { set VTK_DATA "../../../vtkdata" }

# get the interactor ui
source $VTK_TCL/vtkInt.tcl

# create planes
# Create the RenderWindow, Renderer and both Actors
#
vtkRenderer ren1
vtkRenderWindow renWin
    renWin AddRenderer ren1
vtkRenderWindowInteractor iren
    iren SetRenderWindow renWin

# create pipeline
#
vtkPLOT3DReader pl3d
    pl3d SetXYZFileName "$VTK_DATA/combxyz.bin"
    pl3d SetQFileName "$VTK_DATA/combq.bin"
    pl3d SetScalarFunctionNumber 100
    pl3d SetVectorFunctionNumber 202
    pl3d Update

vtkPlaneSource plane
    plane SetResolution 50 50
vtkTransform transP1
    transP1 Translate 3.7 0.0 28.37
    transP1 Scale 5 5 5
    transP1 RotateY 90
vtkTransformPolyDataFilter tpd1
    tpd1 SetInput [plane GetOutput]
    tpd1 SetTransform transP1
vtkOutlineFilter outTpd1
    outTpd1 SetInput [tpd1 GetOutput]
vtkPolyDataMapper mapTpd1
    mapTpd1 SetInput [outTpd1 GetOutput]
vtkActor tpd1Actor
    tpd1Actor SetMapper mapTpd1
    [tpd1Actor GetProperty] SetColor 0 0 0

vtkTransform transP2
    transP2 Translate 9.2 0.0 31.20
    transP2 Scale 5 5 5
    transP2 RotateY 90
vtkTransformPolyDataFilter tpd2
    tpd2 SetInput [plane GetOutput]
    tpd2 SetTransform transP2
vtkOutlineFilter outTpd2
    outTpd2 SetInput [tpd2 GetOutput]
vtkPolyDataMapper mapTpd2
    mapTpd2 SetInput [outTpd2 GetOutput]
vtkActor tpd2Actor
    tpd2Actor SetMapper mapTpd2
    [tpd2Actor GetProperty] SetColor 0 0 0

vtkTransform transP3
    transP3 Translate 13.27 0.0 33.30
    transP3 Scale 5 5 5
    transP3 RotateY 90
vtkTransformPolyDataFilter tpd3
    tpd3 SetInput [plane GetOutput]
    tpd3 SetTransform transP3
vtkOutlineFilter outTpd3
    outTpd3 SetInput [tpd3 GetOutput]
vtkPolyDataMapper mapTpd3
    mapTpd3 SetInput [outTpd3 GetOutput]
vtkActor tpd3Actor
    tpd3Actor SetMapper mapTpd3
    [tpd3Actor GetProperty] SetColor 0 0 0

vtkAppendPolyData appendF
    appendF AddInput [tpd1 GetOutput]
    appendF AddInput [tpd2 GetOutput]
    appendF AddInput [tpd3 GetOutput]

vtkProbeFilter probe
    probe SetInput [appendF GetOutput]
    probe SetSource [pl3d GetOutput]

vtkContourFilter contour
    contour SetInput [probe GetOutput]
    eval contour GenerateValues 50 [[pl3d GetOutput] GetScalarRange]

vtkPolyDataMapper contourMapper
    contourMapper SetInput [contour GetOutput]
    eval contourMapper SetScalarRange [[pl3d GetOutput] GetScalarRange]
vtkActor planeActor
    planeActor SetMapper contourMapper

vtkStructuredGridOutlineFilter outline
    outline SetInput [pl3d GetOutput]
vtkPolyDataMapper outlineMapper
    outlineMapper SetInput [outline GetOutput]
vtkActor outlineActor
    outlineActor SetMapper outlineMapper
    [outlineActor GetProperty] SetColor 0 0 0

ren1 AddActor outlineActor
ren1 AddActor planeActor
ren1 AddActor tpd1Actor
ren1 AddActor tpd2Actor
ren1 AddActor tpd3Actor
ren1 SetBackground 1 1 1
renWin SetSize 500 500

set cam1 [ren1 GetActiveCamera]
$cam1 SetClippingRange 3.95297 50
$cam1 SetFocalPoint 8.88908 0.595038 29.3342
$cam1 SetPosition -12.3332 31.7479 41.2387
$cam1 SetViewUp 0.060772 -0.319905 0.945498
iren Initialize

#renWin SetFileName "probeComb.tcl.ppm"
#renWin SaveImageAsPPM

# render the image
#
iren SetUserMethod {wm deiconify .vtkInteract}

# prevent the tk window from showing up then start the event loop
wm withdraw .



