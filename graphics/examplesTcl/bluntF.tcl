catch {load vtktcl}
# get the interactor ui
source vtkInt.tcl
source "colors.tcl"

# Create the RenderWindow, Renderer and both Actors
#
vtkRenderer ren1
vtkRenderWindow renWin
    renWin AddRenderer ren1
vtkRenderWindowInteractor iren
    iren SetRenderWindow renWin

# read data
#
vtkPLOT3DReader pl3d
    pl3d SetXYZFileName "../../../data/bluntfinxyz.bin"
    pl3d SetQFileName "../../../data/bluntfinq.bin"
    pl3d SetScalarFunctionNumber 100
    pl3d SetVectorFunctionNumber 202
    pl3d Update

# wall
#
vtkStructuredGridGeometryFilter wall
    wall SetInput [pl3d GetOutput]
    wall SetExtent 0 100 0 0 0 100
vtkPolyDataMapper wallMap
    wallMap SetInput [wall GetOutput]
    wallMap ScalarVisibilityOff
vtkActor wallActor
    wallActor SetMapper wallMap
    eval [wallActor GetProperty] SetColor 0.8 0.8 0.8

# fin
# 
vtkStructuredGridGeometryFilter fin
    fin SetInput [pl3d GetOutput]
    fin SetExtent 0 100 0 100 0 0
vtkPolyDataMapper finMap
    finMap SetInput [fin GetOutput]
    finMap ScalarVisibilityOff
vtkActor finActor
    finActor SetMapper finMap
    eval [finActor GetProperty] SetColor 0.8 0.8 0.8

# planes to threshold
vtkStructuredGridGeometryFilter plane1
    plane1 SetInput [pl3d GetOutput]
    plane1 SetExtent 10 10 0 100 0 100
vtkPolyDataMapper plane1Map
    plane1Map SetInput [plane1 GetOutput]
    set pl3dPtData [[pl3d GetOutput] GetPointData]
    set pl3dScalars [$pl3dPtData GetScalars]
    eval plane1Map SetScalarRange [$pl3dScalars GetRange]
vtkActor plane1Actor
    plane1Actor SetMapper plane1Map

vtkStructuredGridGeometryFilter plane2
    plane2 SetInput [pl3d GetOutput]
    plane2 SetExtent 25 25 0 100 0 100
vtkPolyDataMapper plane2Map
    plane2Map SetInput [plane2 GetOutput]
    eval plane2Map SetScalarRange \
      [[[[pl3d GetOutput] GetPointData] GetScalars] GetRange]
vtkActor plane2Actor
    plane2Actor SetMapper plane2Map

vtkStructuredGridGeometryFilter plane3
    plane3 SetInput [pl3d GetOutput]
    plane3 SetExtent 35 35 0 100 0 100
vtkDataSetMapper plane3Map
    plane3Map SetInput [plane3 GetOutput]
    eval plane3Map SetScalarRange \
      [[[[pl3d GetOutput] GetPointData] GetScalars] GetRange]
vtkActor plane3Actor
    plane3Actor SetMapper plane3Map

# outline
vtkStructuredGridOutlineFilter outline
    outline SetInput [pl3d GetOutput]
vtkPolyDataMapper outlineMapper
    outlineMapper SetInput [outline GetOutput]
vtkActor outlineActor
    outlineActor SetMapper outlineMapper
    set outlineProp [outlineActor GetProperty]
    eval $outlineProp SetColor 0 0 0

# Add the actors to the renderer, set the background and size
#
ren1 AddActor outlineActor
ren1 AddActor wallActor
ren1 AddActor finActor
ren1 AddActor plane1Actor
ren1 AddActor plane2Actor
ren1 AddActor plane3Actor
ren1 SetBackground 1 1 1
renWin SetSize 500 500

set cam1 [ren1 GetActiveCamera]
$cam1 Azimuth -40
$cam1 Zoom 1.4

iren Initialize
renWin Render

# render the image
#
iren SetUserMethod {wm deiconify .vtkInteract}

# prevent the tk window from showing up then start the event loop
wm withdraw .

#renWin SetFileName bluntF.tcl.ppm
#renWin SaveImageAsPPM
