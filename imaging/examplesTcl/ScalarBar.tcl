# This example demonstrates how to use cell data as well as the programmable attribute 
# filter. Example randomly colors cells with scalar values.
catch {load vtktcl}
if { [catch {set VTK_TCL $env(VTK_TCL)}] != 0} { set VTK_TCL "../../examplesTcl" }
if { [catch {set VTK_DATA $env(VTK_DATA)}] != 0} { set VTK_DATA "../../../vtkdata" }


# get the interactor ui
source $VTK_TCL/vtkInt.tcl

# create pipeline
#
# create sphere to color
vtkSphereSource sphere
    sphere SetThetaResolution 20
    sphere SetPhiResolution 40

# Compute random scalars (colors) for each cell
vtkProgrammableAttributeDataFilter randomColors
    randomColors SetInput [sphere GetOutput]
    randomColors SetExecuteMethod colorCells

proc colorCells {} {
    vtkMath randomColorGenerator
    set input [randomColors GetInput]
    set output [randomColors GetOutput]
    set numCells [$input GetNumberOfCells]
    vtkScalars colors
	colors SetNumberOfScalars $numCells

    for {set i 0} {$i < $numCells} {incr i} {
        colors SetScalar $i [randomColorGenerator Random 0 1]
    }

    [$output GetCellData] CopyScalarsOff
    [$output GetCellData] PassData [$input GetCellData]
    [$output GetCellData] SetScalars colors

    colors Delete; #reference counting - it's ok
    randomColorGenerator Delete
}
# mapper and actor
vtkPolyDataMapper mapper
    mapper SetInput [randomColors GetPolyDataOutput]
    eval mapper SetScalarRange [[randomColors GetPolyDataOutput] GetScalarRange]
vtkActor sphereActor
    sphereActor SetMapper mapper

# Create a scalar bar
vtkScalarBarActor scalarBar
    scalarBar SetLookupTable [mapper GetLookupTable]
    scalarBar SetTitle "Temperature"
    [scalarBar GetPositionCoordinate] SetCoordinateSystemToNormalizedViewport
    [scalarBar GetPositionCoordinate] SetValue 0.1 0.01
    scalarBar SetOrientationToHorizontal
    scalarBar SetWidth 0.8
    scalarBar SetHeight 0.17

# Create graphics stuff
# Create the RenderWindow, Renderer and both Actors
#
vtkRenderer ren1
vtkRenderWindow renWin
    renWin AddRenderer ren1
vtkRenderWindowInteractor iren
    iren SetRenderWindow renWin

ren1 AddActor sphereActor
ren1 AddActor2D scalarBar
renWin SetSize 500 500

# render the image
#
iren SetUserMethod {wm deiconify .vtkInteract}
[ren1 GetActiveCamera] Zoom 1.5
renWin Render
scalarBar SetNumberOfLabels 8
renWin Render

#renWin SetFileName "ScalarBar.tcl.ppm"
#renWin SaveImageAsPPM

# prevent the tk window from showing up then start the event loop
wm withdraw .

