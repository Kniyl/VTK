package require vtk
package require vtkinteraction

# Read some AVS UCD data in ASCII form
vtkAVSucdReader r
   r SetFileName "$VTK_DATA_ROOT/Data/AVS/cellsnd.ascii.inp"
vtkDataSetMapper AVSMapper
   AVSMapper SetInput [r GetOutput]
vtkActor AVSActor
   AVSActor SetMapper AVSMapper

# Read some AVS UCD data in binary form
vtkAVSucdReader r2
   r2 SetFileName "$VTK_DATA_ROOT/Data/AVS/cellsnd.bin.inp"
vtkDataSetMapper AVSMapper2
   AVSMapper2 SetInput [r2 GetOutput]
vtkActor AVSActor2
   AVSActor2 SetMapper AVSMapper2
   AVSActor2 AddPosition 10 0 0

# Create the RenderWindow, Renderer and both Actors
#
vtkRenderer ren1
vtkRenderWindow renWin
    renWin AddRenderer ren1
vtkRenderWindowInteractor iren
    iren SetRenderWindow renWin

# Add the actors to the renderer, set the background and size
#
ren1 AddActor AVSActor
ren1 AddActor AVSActor2

renWin SetSize 300 150
iren Initialize
renWin Render
[ren1 GetActiveCamera] Zoom 3.0

iren AddObserver UserEvent {wm deiconify .vtkInteract}

# prevent the tk window from showing up then start the event loop
wm withdraw .
