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

# create spring profile (a circle)
#
vtkFloatPoints points
    points InsertPoint 0 1.0 0.0 0.0
    points InsertPoint 1 1.0732 0.0 -0.1768
    points InsertPoint 2 1.25 0.0 -0.25
    points InsertPoint 3 1.4268 0.0 -0.1768
    points InsertPoint 4 1.5 0.0 0.00
    points InsertPoint 5 1.4268 0.0 0.1768
    points InsertPoint 6 1.25 0.0 0.25
    points InsertPoint 7 1.0732 0.0 0.1768
vtkCellArray poly
    poly InsertNextCell 8;#number of points
    poly InsertCellPoint 0
    poly InsertCellPoint 1
    poly InsertCellPoint 2
    poly InsertCellPoint 3
    poly InsertCellPoint 4
    poly InsertCellPoint 5
    poly InsertCellPoint 6
    poly InsertCellPoint 7
vtkPolyData profile
    profile SetPoints points
    profile SetPolys poly

# extrude profile to make spring
#
vtkRotationalExtrusionFilter extrude
    extrude SetInput profile
    extrude SetResolution 360
    extrude SetTranslation 6
    extrude SetDeltaRadius 1.0
    extrude SetAngle 2160.0;#six revolutions
    
vtkPolyDataNormals normals
    normals SetInput [extrude GetOutput]
    normals SetFeatureAngle 60

vtkPolyDataMapper map
    map SetInput [normals GetOutput]

vtkActor spring
    spring SetMapper map
    [spring GetProperty] SetColor 0.6902 0.7686 0.8706
    [spring GetProperty] SetDiffuse 0.7
    [spring GetProperty] SetSpecular 0.4
    [spring GetProperty] SetSpecularPower 20
    [spring GetProperty] BackfaceCullingOn

# Add the actors to the renderer, set the background and size
#
ren1 AddActor spring
ren1 SetBackground 1 1 1
renWin SetSize 500 500

set cam1 [ren1 GetActiveCamera]
$cam1 Azimuth 90

renWin Render
#renWin SetFileName "spring.tcl.ppm"
#renWin SaveImageAsPPM

iren SetUserMethod {wm deiconify .vtkInteract}

# prevent the tk window from showing up then start the event loop
wm withdraw .



