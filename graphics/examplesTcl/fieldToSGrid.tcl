catch {load vtktcl}
if { [catch {set VTK_TCL $env(VTK_TCL)}] != 0} { set VTK_TCL "../../examplesTcl" }
if { [catch {set VTK_DATA $env(VTK_DATA)}] != 0} { set VTK_DATA "../../../vtkdata" }


# Demonstrate the generation of a structured grid from field data. The output
# should be similar to combIso.tcl.

# get the interactor ui
source $VTK_TCL/vtkInt.tcl
source $VTK_TCL/colors.tcl

# Create a reader and write out the field
vtkPLOT3DReader pl3d
    pl3d SetXYZFileName "$VTK_DATA/combxyz.bin"
    pl3d SetQFileName "$VTK_DATA/combq.bin"
vtkDataSetToDataObjectFilter ds2do
    ds2do SetInput [pl3d GetOutput]
vtkDataObjectWriter writer
    writer SetInput [ds2do GetOutput]
    writer SetFileName "SGridField.vtk"
    writer Write

# read the field
vtkDataObjectReader dor
    dor SetFileName "SGridField.vtk"
vtkDataObjectToDataSetFilter do2ds
    do2ds SetInput [dor GetOutput]
    do2ds SetDataSetTypeToStructuredGrid
    do2ds SetDimensionsComponent Dimensions 0 
    do2ds SetPointComponent 0 Points 0 
    do2ds SetPointComponent 1 Points 1 
    do2ds SetPointComponent 2 Points 2 
vtkFieldDataToAttributeDataFilter fd2ad
    fd2ad SetInput [do2ds GetStructuredGridOutput]
    fd2ad SetInputFieldToDataObjectField
    fd2ad SetOutputAttributeDataToPointData
    fd2ad SetVectorComponent 0 PointVectors 0 
    fd2ad SetVectorComponent 1 PointVectors 1 
    fd2ad SetVectorComponent 2 PointVectors 2 
    fd2ad SetScalarComponent 0 PointScalars 0 

# create pipeline
#
vtkContourFilter iso
    iso SetInput [fd2ad GetOutput]
    iso SetValue 0 .38
vtkPolyDataNormals normals
    normals SetInput [iso GetOutput]
    normals SetFeatureAngle 45
vtkPolyDataMapper isoMapper
    isoMapper SetInput [normals GetOutput]
    isoMapper ScalarVisibilityOff
vtkActor isoActor
    isoActor SetMapper isoMapper
    eval [isoActor GetProperty] SetColor $bisque

vtkStructuredGridOutlineFilter outline
    outline SetInput [fd2ad GetStructuredGridOutput]
vtkPolyDataMapper outlineMapper
    outlineMapper SetInput [outline GetOutput]
vtkActor outlineActor
    outlineActor SetMapper outlineMapper

# Create the RenderWindow, Renderer and both Actors
#
vtkRenderer ren1
vtkRenderWindow renWin
    renWin AddRenderer ren1
vtkRenderWindowInteractor iren
    iren SetRenderWindow renWin

# Add the actors to the renderer, set the background and size
#
ren1 AddActor outlineActor
ren1 AddActor isoActor
ren1 SetBackground 1 1 1
renWin SetSize 500 500
ren1 SetBackground 0.1 0.2 0.4

set cam1 [ren1 GetActiveCamera]
$cam1 SetClippingRange 3.95297 50
$cam1 SetFocalPoint 9.71821 0.458166 29.3999
$cam1 SetPosition 2.7439 -37.3196 38.7167
$cam1 SetViewUp -0.16123 0.264271 0.950876

# render the image
#
iren SetUserMethod {wm deiconify .vtkInteract}

renWin Render
renWin SetFileName "fieldToSGrid.tcl.ppm"
#renWin SaveImageAsPPM

# prevent the tk window from showing up then start the event loop
wm withdraw .


