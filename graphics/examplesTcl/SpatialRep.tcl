catch {load vtktcl}
# get the interactor ui
source ../../examplesTcl/vtkInt.tcl

# lines make a nice test
vtkLineSource line1
  line1 SetPoint1 0 0 0
  line1 SetPoint2 1 0 0
  line1 SetResolution 1000
vtkLineSource line2
  line2 SetPoint1 0 0 0
  line2 SetPoint2 1 1 1
  line2 SetResolution 1000
#vtkAppendPolyData asource
#  asource AddInput [line1 GetOutput]
#  asource AddInput [line2 GetOutput]

vtkSTLReader asource
  asource SetFileName ../../../vtkdata/42400-IDGH.stl
#vtkCyberReader asource
#  asource SetFileName ../../../vtkdata/fran_cut
vtkPolyDataMapper dataMapper
  dataMapper SetInput [asource GetOutput]
vtkActor model
  model SetMapper dataMapper
  [model GetProperty] SetColor 1 0 0
#  model VisibilityOff

#vtkPointLocator locator
#vtkOBBTree locator
vtkCellLocator locator
  locator SetMaxLevel 4
  locator AutomaticOff
vtkSpatialRepresentationFilter boxes
  boxes SetInput [asource GetOutput]
  boxes SetSpatialRepresentation locator
vtkPolyDataMapper boxMapper
  boxMapper SetInput [boxes GetOutput]
#  boxMapper SetInput [boxes GetOutput 2]
vtkActor boxActor
  boxActor SetMapper boxMapper
  [boxActor GetProperty] SetDiffuse 0
  [boxActor GetProperty] SetAmbient 1
  [boxActor GetProperty] SetRepresentationToWireframe

vtkRenderer ren1
vtkRenderWindow renWin
    renWin AddRenderer ren1
vtkRenderWindowInteractor iren
    iren SetRenderWindow renWin

# Add the actors to the renderer, set the background and size
#
ren1 AddActor model
ren1 AddActor boxActor
ren1 SetBackground 0.1 0.2 0.4
renWin SetSize 500 500
iren Initialize

# render the image
#
iren SetUserMethod {wm deiconify .vtkInteract}
[ren1 GetActiveCamera] Zoom 1.4
renWin Render

#renWin SetFileName valid/SpatialRep.tcl.ppm
#renWin SaveImageAsPPM

# prevent the tk window from showing up then start the event loop
wm withdraw .

