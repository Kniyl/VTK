catch {load vtktcl}
# get the interactor ui
source vtkInt.tcl
source "colors.tcl"
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
    pl3d SetXYZFileName "../../../data/combxyz.bin"
    pl3d SetQFileName "../../../data/combq.bin"
    pl3d SetScalarFunctionNumber 100
    pl3d SetVectorFunctionNumber 202
    pl3d Update
vtkStructuredGridGeometryFilter plane
    plane SetInput [pl3d GetOutput]
    plane SetExtent 0 100 0 100 0 0
vtkPolyDataMapper planeMapper
    planeMapper SetInput [plane GetOutput]
    planeMapper SetScalarRange 0.197813 0.710419
vtkActor planeActor
    planeActor SetMapper planeMapper

vtkStructuredGridOutlineFilter outline
    outline SetInput [pl3d GetOutput]
vtkPolyDataMapper outlineMapper
    outlineMapper SetInput [outline GetOutput]
vtkActor outlineActor
    outlineActor SetMapper outlineMapper
    eval [outlineActor GetProperty] SetColor $black

# Add the actors to the renderer, set the background and size
#
ren1 AddActor outlineActor
ren1 AddActor planeActor
ren1 SetBackground 1 1 1
renWin SetSize 500 500
iren Initialize

set cam1 [ren1 GetActiveCamera]
    $cam1 SetClippingRange 3.95297 50
    $cam1 SetFocalPoint 8.88908 0.595038 29.3342
    $cam1 SetPosition -12.3332 31.7479 41.2387
    $cam1 ComputeViewPlaneNormal
    $cam1 SetViewUp 0.060772 -0.319905 0.945498

# render the image
#
iren SetUserMethod {wm deiconify .vtkInteract}

for {set j 0} {$j<3} {incr j 1} {
    for {set i 0} {$i<25} {incr i 1} {
        eval plane SetExtent 0 100 0 100 $i $i
        renWin Render
    }
}

#renWin SetFileName "movePlane.tcl.ppm"
#renWin SaveImageAsPPM

# prevent the tk window from showing up then start the event loop
wm withdraw .



