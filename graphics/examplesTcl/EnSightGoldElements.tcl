catch {load vtktcl}
if { [catch {set VTK_TCL $env(VTK_TCL)}] != 0} { set VTK_TCL "../../examplesTcl" }
if { [catch {set VTK_DATA $env(VTK_DATA)}] != 0} { set VTK_DATA "../../../vtkdata" }

# user interface command widget
source $VTK_TCL/vtkInt.tcl

# create a rendering window and renderer
vtkRenderer ren1
vtkRenderWindow renWin
    renWin AddRenderer ren1
    renWin StereoCapableWindowOn  
vtkRenderWindowInteractor iren
    iren SetRenderWindow renWin

vtkGenericEnSightReader reader
    reader SetCaseFileName "C:/EnSightFiles/elements/elements.case"
    reader Update

vtkGeometryFilter geom0
    geom0 SetInput [reader GetOutput]
vtkGeometryFilter geom1
    geom1 SetInput [reader GetOutput 1]

vtkPolyDataMapper mapper0
    mapper0 SetInput [geom0 GetOutput]
    mapper0 SetColorModeToMapScalars
    mapper0 SetScalarModeToUsePointFieldData
    mapper0 ColorByArrayComponent pointScalars 0
    mapper0 SetScalarRange 0 112
vtkPolyDataMapper mapper1
    mapper1 SetInput [geom1 GetOutput]
    mapper1 SetColorModeToMapScalars
    mapper1 SetScalarModeToUsePointFieldData
    mapper1 ColorByArrayComponent pointScalars 0
    mapper1 SetScalarRange 0 112

vtkActor actor0
    actor0 SetMapper mapper0
vtkActor actor1
    actor1 SetMapper mapper1

# assign our actor to the renderer
ren1 AddActor actor0
ren1 AddActor actor1

# enable user interface interactor
iren SetUserMethod {wm deiconify .vtkInteract}
iren Initialize

# prevent the tk window from showing up then start the event loop
wm withdraw .

