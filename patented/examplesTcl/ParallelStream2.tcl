# This example test the ParallelStreaming flag in the 
# vtkAppendPolyData filter.

# parameters
set NUMBER_OF_PIECES 7



catch {load vtktcl}
# get the interactor ui
source ../../examplesTcl/vtkInt.tcl
source ../../examplesTcl/colors.tcl

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
    reader SetFilePrefix "../../../vtkdata/headsq/half"
    reader SetDataSpacing 1.6 1.6 1.5

vtkAppendPolyData app
    app ParallelStreamingOn

vtkAppendPolyData appa
    appa ParallelStreamingOn

vtkAppendPolyData appb
    appb ParallelStreamingOn


for {set i 0} {$i < $NUMBER_OF_PIECES} {incr i} {

  vtkSynchronizedTemplates3D iso$i
    iso$i SetInput [reader GetOutput]
    iso$i SetValue 0 500
    iso$i ComputeScalarsOff
    iso$i ComputeGradientsOff

  set val [expr 0.0 + $i / ($NUMBER_OF_PIECES - 1.0)]
  vtkElevationFilter elev$i
    elev$i SetInput [iso$i GetOutput]
    elev$i SetScalarRange $val [expr $val+0.001]

  appa AddInput [elev$i GetOutput]
  appb AddInput [elev$i GetOutput]
}

app AddInput [appa GetOutput]
app AddInput [appb GetOutput]


vtkPolyDataMapper mapper
  mapper SetInput [app GetOutput]
  mapper ImmediateModeRenderingOn

vtkActor actor
  actor SetMapper mapper
  [actor GetProperty] SetSpecularPower 30
  [actor GetProperty] SetDiffuse .7
  [actor GetProperty] SetSpecular .5
    
ren1 AddActor actor

vtkOutlineFilter outline
  outline SetInput [reader GetOutput]
vtkPolyDataMapper outlineMapper
  outlineMapper SetInput [outline GetOutput]
vtkActor outlineActor
  outlineActor SetMapper outlineMapper
  outlineActor VisibilityOff

# Add the actors to the renderer, set the background and size
#
ren1 AddActor outlineActor
ren1 SetBackground 0.9 .9 .9
[ren1 GetActiveCamera] SetFocalPoint 100 100 65
[ren1 GetActiveCamera] SetPosition 100 450 65
[ren1 GetActiveCamera] SetViewUp 0 0 -1
[ren1 GetActiveCamera] SetViewAngle 30
[ren1 GetActiveCamera] ComputeViewPlaneNormal
ren1 ResetCameraClippingRange


renWin SetSize 450 450
iren Initialize

# render the image
#
iren SetUserMethod {wm deiconify .vtkInteract}
renWin SetFileName "ParallelStream2.tcl.ppm"
#renWin SaveImageAsPPM

# prevent the tk window from showing up then start the event loop
wm withdraw .


