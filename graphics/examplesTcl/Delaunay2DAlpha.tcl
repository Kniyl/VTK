catch {load vtktcl}
if { [catch {set VTK_TCL $env(VTK_TCL)}] != 0} { set VTK_TCL "../../examplesTcl" }
if { [catch {set VTK_DATA $env(VTK_DATA)}] != 0} { set VTK_DATA "../../../vtkdata" }

# get the interactor ui
source $VTK_TCL/vtkInt.tcl

# Create the RenderWindow, Renderer and both Actors
#
vtkRenderer ren1
vtkRenderWindow renWin
    renWin AddRenderer ren1
vtkRenderWindowInteractor iren
    iren SetRenderWindow renWin

# create some points
#
vtkMath math
vtkPoints points
for {set i 0} {$i<1000} {incr i 1} {
    eval points InsertPoint $i [math Random 0 1] [math Random 0 1] 0.0
}

vtkPolyData profile
    profile SetPoints points

# triangulate them
#
vtkDelaunay2D del
    del SetInput profile
    del BoundingTriangulationOn
    del SetTolerance 0.001
    del SetAlpha .1
    del Update
    
vtkShrinkPolyData shrink
    shrink SetInput [del GetOutput]

vtkPolyDataMapper map
    map SetInput [shrink GetOutput]

vtkActor triangulation
    triangulation SetMapper map
    [triangulation GetProperty] SetColor 1 0 0

# Add the actors to the renderer, set the background and size
#
ren1 AddActor triangulation
ren1 SetBackground 1 1 1
renWin SetSize 500 500
renWin Render

set cam1 [ren1 GetActiveCamera]
$cam1 Zoom 1.5

# render the image
#
iren SetUserMethod {wm deiconify .vtkInteract}

renWin Render

#renWin SetFileName Delaunay2DAlpha.tcl.ppm
#renWin SaveImageAsPPM

# prevent the tk window from showing up then start the event loop
wm withdraw .
