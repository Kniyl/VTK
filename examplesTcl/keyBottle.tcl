#
# Show keyframing of camera
#
catch {load vtktcl}
# get the interactor ui
source vtkInt.tcl

# Create the RenderWindow, Renderer and both Actors
#
vtkRenderer ren1
vtkRenderWindow renWin
    renWin AddRenderer ren1
vtkRenderWindowInteractor iren
    iren SetRenderWindow renWin

# create bottle profile
#
vtkFloatPoints points
  points InsertPoint 0 0.01 0.0 0.0
  points InsertPoint 1 1.5 0.0 0.0
  points InsertPoint 2 1.5 0.0 3.5
  points InsertPoint 3 1.25 0.0 3.75
  points InsertPoint 4 0.75 0.0 4.00
  points InsertPoint 5 0.6 0.0 4.35
  points InsertPoint 6 0.7 0.0 4.65
  points InsertPoint 7 1.0 0.0 4.75
  points InsertPoint 8 1.0 0.0 5.0
  points InsertPoint 9 0.01 0.0 5.0
vtkCellArray lines
  lines InsertNextCell 10;#number of points
  lines InsertCellPoint 0
  lines InsertCellPoint 1
  lines InsertCellPoint 2
  lines InsertCellPoint 3
  lines InsertCellPoint 4
  lines InsertCellPoint 5
  lines InsertCellPoint 6
  lines InsertCellPoint 7
  lines InsertCellPoint 8
  lines InsertCellPoint 9
vtkPolyData profile
    profile SetPoints points
    profile SetLines lines

# extrude profile to make bottle
#
vtkRotationalExtrusionFilter extrude
  extrude SetInput profile
  extrude SetResolution 60
    
vtkPolyDataMapper map
    map SetInput [extrude GetOutput]

vtkActor bottle
    bottle SetMapper map
    [bottle GetProperty] SetColor 0.3800 0.7000 0.1600

# Add the actors to the renderer, set the background and size
#
ren1 AddActor bottle
ren1 SetBackground 1 1 1

renWin SetSize 500 500
renWin Render

set cam1 [ren1 GetActiveCamera]
$cam1 Elevation 30

# render the image
#
iren SetUserMethod {wm deiconify .vtkInteract}
iren Initialize

# prevent the tk window from showing up then start the event loop
wm withdraw .

#
# get keyframe procs
source KeyFrame.tcl
#
# keyframe the azimuth
source azimuthKey.tcl



