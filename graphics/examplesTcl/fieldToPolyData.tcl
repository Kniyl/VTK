catch {load vtktcl}

# This example demonstrates the reading of a field and conversion to PolyData
# The output should be the same as polyEx.tcl.

# get the interactor ui
source ../../examplesTcl/vtkInt.tcl

# Create a reader and write out the field
vtkPolyDataReader reader
    reader SetFileName "../../../vtkdata/polyEx.vtk"
vtkDataSetToDataObjectFilter ds2do
    ds2do SetInput [reader GetOutput]
vtkDataObjectWriter writer
    writer SetInput [ds2do GetOutput]
    writer SetFileName "PolyField.vtk"
    writer Write

# create pipeline
#
vtkDataObjectReader dor
    dor SetFileName "PolyField.vtk"
vtkDataObjectToDataSetFilter do2ds
    do2ds SetInput [dor GetOutput]
    do2ds SetDataSetTypeToPolyData
    do2ds SetPointComponent 0 Points 0 
    do2ds SetPointComponent 1 Points 1 
    do2ds SetPointComponent 2 Points 2 
    do2ds SetPolysComponent Polys 0
vtkFieldDataToAttributeDataFilter fd2ad
    fd2ad SetInput [do2ds GetPolyDataOutput]
    fd2ad SetInputFieldToDataObjectField
    fd2ad SetOutputAttributeDataToPointData
    fd2ad SetScalarComponent 0 PointScalars 0 

vtkPolyDataMapper mapper
    mapper SetInput [fd2ad GetPolyDataOutput]
    eval mapper SetScalarRange [[fd2ad GetOutput] GetScalarRange]
vtkActor actor
    actor SetMapper mapper

# Create the RenderWindow, Renderer and both Actors
vtkRenderer ren1
vtkRenderWindow renWin
    renWin AddRenderer ren1
vtkRenderWindowInteractor iren
    iren SetRenderWindow renWin

ren1 AddActor actor
ren1 SetBackground 1 1 1
renWin SetSize 500 500

set cam1 [ren1 GetActiveCamera]
$cam1 SetClippingRange .348 17.43
$cam1 SetPosition 2.92 2.62 -0.836
$cam1 ComputeViewPlaneNormal
$cam1 SetViewUp -0.436 -0.067 -0.897

# render the image
#
iren SetUserMethod {wm deiconify .vtkInteract}
renWin Render

renWin SetFileName "fieldToPolyData.tcl.ppm"
#renWin SaveImageAsPPM

# prevent the tk window from showing up then start the event loop
wm withdraw .


