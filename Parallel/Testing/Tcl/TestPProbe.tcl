# ParaView Version 0.5

package require vtk

if { [ info command vtkMesaRenderer ] != "" } {
    vtkGraphicsFactory _graphics_fact
    _graphics_fact SetUseMesaClasses 1
    _graphics_fact Delete
}

# create a rendering window and renderer
vtkRenderer Ren1
	Ren1 SetBackground .5 .8 1
vtkRenderWindow renWin
	renWin AddRenderer Ren1
	renWin SetSize 529 586
        if { $myProcId > 0 } {
           renWin OffScreenRenderingOn
        }

# camera parameters
set camera [Ren1 GetActiveCamera]
	$camera SetPosition 199.431 196.879 15.7781
	$camera SetFocalPoint 33.5 33.5 33.5
	$camera SetViewUp 0.703325 -0.702557 0.108384
	$camera SetViewAngle 30
	$camera SetClippingRange 132.14 361.741

vtkPDataSetReader ironProt0
	ironProt0 SetFileName "$VTK_DATA_ROOT/Data/ironProt.vtk"

vtkPVGeometryFilter Geometry4
	Geometry4 SetInput [ironProt0 GetOutput]
vtkPolyDataMapper Mapper4
	Mapper4 SetInput [Geometry4 GetOutput]
	Mapper4 SetImmediateModeRendering 0
	Mapper4 SetScalarRange 0 1
	Mapper4 SetScalarVisibility 0
	Mapper4 SetScalarModeToDefault
vtkActor Actor4
	Actor4 SetMapper Mapper4
	[Actor4 GetProperty] SetRepresentationToSurface
	[Actor4 GetProperty] SetInterpolationToGouraud
[Actor4 GetProperty] SetColor 1 1 1
Ren1 AddActor Actor4

vtkLineSource probeLine
  probeLine SetPoint1 0 67 10
  probeLine SetPoint2 67 0 50
  probeLine SetResolution 500

vtkMultiProcessController controler

vtkPProbeFilter Probe0
	Probe0 SetSource [ironProt0 GetOutput]
        Probe0 SetInput [ probeLine GetOutput ]
        Probe0 SetController controler

vtkTubeFilter Tuber0
 	Tuber0 SetInput [Probe0 GetOutput]
 	Tuber0 SetNumberOfSides 10
 	Tuber0 SetCapping 0
 	Tuber0 SetRadius 1
 	Tuber0 SetVaryRadius 1
 	Tuber0 SetRadiusFactor 10
vtkPVGeometryFilter Geometry6
 	Geometry6 SetInput [Tuber0 GetOutput]
Geometry6 Update

vtkPolyDataMapper Mapper6
 	Mapper6 SetInput [Geometry6 GetOutput]
 	Mapper6 SetImmediateModeRendering 0
 	Mapper6 SetScalarRange 0 228
 	Mapper6 SetScalarVisibility 1
 	Mapper6 SetScalarModeToUsePointFieldData
 	Mapper6 ColorByArrayComponent {scalars} -1
        Mapper6 UseLookupTableScalarRangeOn
vtkActor Actor6
 	Actor6 SetMapper Mapper6
 	[Actor6 GetProperty] SetRepresentationToSurface
 	[Actor6 GetProperty] SetInterpolationToGouraud
Ren1 AddActor Actor6

if { $numProcs > 1 } {
    compManager SetRenderWindow renWin 
        compManager InitializePieces
}

renWin SetWindowName "Process $myProcId"

wm withdraw .

renWin Render

