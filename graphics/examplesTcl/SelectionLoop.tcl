catch {load vtktcl}
if { [catch {set VTK_TCL $env(VTK_TCL)}] != 0} { set VTK_TCL "../../examplesTcl" }
if { [catch {set VTK_DATA $env(VTK_DATA)}] != 0} { set VTK_DATA "../../../vtkdata" }

#
# Demonstrate the use of implicit selection loop as well as closest point
# connectivity
#
source $VTK_TCL/vtkInt.tcl
source $VTK_TCL/colors.tcl

# create pipeline
#
vtkSphereSource sphere
    sphere SetRadius 1
    sphere SetPhiResolution 100
    sphere SetThetaResolution 100
vtkPoints selectionPoints
    selectionPoints InsertPoint 0 0.07325 0.8417 0.5612
    selectionPoints InsertPoint 1 0.07244 0.6568 0.7450
    selectionPoints InsertPoint 2 0.1727 0.4597 0.8850
    selectionPoints InsertPoint 3 0.3265 0.6054 0.7309
    selectionPoints InsertPoint 4 0.5722 0.5848 0.5927
    selectionPoints InsertPoint 5 0.4305 0.8138 0.4189
vtkImplicitSelectionLoop loop
    loop SetLoop selectionPoints
vtkExtractGeometry extract
    extract SetInput [sphere GetOutput]
    extract SetImplicitFunction loop
vtkConnectivityFilter connect
    connect SetInput [extract GetOutput]
    connect SetExtractionModeToClosestPointRegion
    eval connect SetClosestPoint [selectionPoints GetPoint 0]
vtkDataSetMapper clipMapper
    clipMapper SetInput [connect GetOutput]
vtkProperty backProp
    eval backProp SetDiffuseColor $tomato
vtkActor clipActor
    clipActor SetMapper clipMapper
    eval [clipActor GetProperty] SetColor $peacock
    clipActor SetBackfaceProperty backProp

# Create graphics stuff
#
vtkRenderer ren1
vtkRenderWindow renWin
    renWin AddRenderer ren1
vtkRenderWindowInteractor iren
    iren SetRenderWindow renWin

# Add the actors to the renderer, set the background and size
#
ren1 AddActor clipActor
ren1 SetBackground 1 1 1
[ren1 GetActiveCamera] Azimuth 30
[ren1 GetActiveCamera] Elevation 30
[ren1 GetActiveCamera] Dolly 1.2

renWin SetSize 400 400
iren Initialize

# render the image
#
iren SetUserMethod {wm deiconify .vtkInteract}
renWin SetFileName "SelectionLoop.tcl.ppm"
#renWin SaveImageAsPPM

# prevent the tk window from showing up then start the event loop
wm withdraw .


