catch {load vtktcl}
#create triangle strip - won't see anything with backface culling on

# get the interactor ui
source ../../examplesTcl/vtkInt.tcl

# Create the RenderWindow, Renderer and both Actors
#
vtkRenderer ren1
vtkRenderWindow renWin
    renWin AddRenderer ren1
vtkRenderWindowInteractor iren
    iren SetRenderWindow renWin

# create triangle strip
#
vtkPoints points
    points InsertPoint 0 0.0 0.0 0.0
    points InsertPoint 1 0.0 1.0 0.0
    points InsertPoint 2 1.0 0.0 0.0
    points InsertPoint 3 1.0 1.0 0.0
    points InsertPoint 4 2.0 0.0 0.0
    points InsertPoint 5 2.0 1.0 0.0
    points InsertPoint 6 3.0 0.0 0.0
    points InsertPoint 7 3.0 1.0 0.0
vtkCellArray strips
    strips InsertNextCell 8;#number of points
    strips InsertCellPoint 0
    strips InsertCellPoint 1
    strips InsertCellPoint 2
    strips InsertCellPoint 3
    strips InsertCellPoint 4
    strips InsertCellPoint 5
    strips InsertCellPoint 6
    strips InsertCellPoint 7
vtkPolyData profile
    profile SetPoints points
    profile SetStrips strips

vtkPolyDataMapper map
    map SetInput profile

vtkActor strip
    strip SetMapper map
    [strip GetProperty] SetColor 0.3800 0.7000 0.1600
    [strip GetProperty] BackfaceCullingOff

# Add the actors to the renderer, set the background and size
#
ren1 AddActor strip
ren1 SetBackground 1 1 1
renWin SetSize 500 500
renWin Render

#renWin SetFileName "strip.tcl.ppm"
#renWin SaveImageAsPPM

# render the image
#
iren SetUserMethod {wm deiconify .vtkInteract}

# prevent the tk window from showing up then start the event loop
wm withdraw .



