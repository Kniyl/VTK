package require vtktcl_interactor

vtkParallelFactory pf
pf RegisterFactory pf

vtkSphereSource sphere
sphere SetPhiResolution 32
sphere SetThetaResolution 32

vtkExtractPolyDataPiece extract
extract SetInput [sphere GetOutput]

vtkPolyDataNormals normals
normals SetInput [extract GetOutput]

vtkPieceScalars ps
ps SetInput [normals GetOutput]

vtkPolyDataMapper mapper
mapper SetInput [ps GetOutput]
mapper SetNumberOfPieces 2

vtkActor actor
actor SetMapper mapper

vtkSphereSource sphere2
sphere2 SetPhiResolution 32
sphere2 SetThetaResolution 32

vtkExtractPolyDataPiece extract2
extract2 SetInput [sphere2 GetOutput]

vtkPolyDataMapper mapper2
mapper2 SetInput [extract2 GetOutput]
mapper2 SetNumberOfPieces 2
mapper2 SetPiece 1
mapper2 SetScalarRange 0 4
mapper2 SetScalarModeToUseCellFieldData
mapper2 SetColorModeToMapScalars 
mapper2 ColorByArrayComponent "vtkGhostLevels" 0
mapper2 SetGhostLevel 4

vtkActor actor2
actor2 SetMapper mapper2
actor2 SetPosition 1.5 0 0

vtkSphereSource sphere3
sphere3 SetPhiResolution 32
sphere3 SetThetaResolution 32

vtkExtractPolyDataPiece extract3
extract3 SetInput [sphere3 GetOutput]

vtkPieceScalars ps3
ps3 SetInput [extract3 GetOutput]

vtkPolyDataMapper mapper3
mapper3 SetInput [ps3 GetOutput]
mapper3 SetNumberOfSubPieces 8
mapper3 SetScalarRange 0 8

vtkActor actor3
actor3 SetMapper mapper3
actor3 SetPosition 0 -1.5 0

vtkRenderer ren
ren AddActor actor
ren AddActor actor2
ren AddActor actor3

vtkRenderWindow renWin
renWin AddRenderer ren

vtkRenderWindowInteractor iren
iren SetRenderWindow renWin
iren Initialize
iren SetUserMethod {wm deiconify .vtkInteract}

wm withdraw .

