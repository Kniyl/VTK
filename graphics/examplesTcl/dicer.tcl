catch {load vtktcl}
# get the interactor ui and colors
source vtkInt.tcl
source "colors.tcl"

# Create the RenderWindow, Renderer and both Actors
#
vtkRenderer ren1
vtkRenderWindow renWin
    renWin AddRenderer ren1
vtkRenderWindowInteractor iren
    iren SetRenderWindow renWin

# create pipeline
#
vtkSTLReader reader
    reader SetFileName ../../../data/42400-IDGH.stl
vtkDicer dicer
    dicer SetInput [reader GetOutput]
    dicer SetNumberOfPointsPerPiece 1000
    dicer Update
vtkDataSetMapper isoMapper
    isoMapper SetInput [dicer GetOutput]
    isoMapper SetScalarRange 0 [dicer GetNumberOfPieces]
vtkActor isoActor
    isoActor SetMapper isoMapper
    eval [isoActor GetProperty] SetColor $raw_sienna

vtkOutlineFilter outline
    outline SetInput [reader GetOutput]
vtkPolyMapper outlineMapper
    outlineMapper SetInput [outline GetOutput]
vtkActor outlineActor
    outlineActor SetMapper outlineMapper
    [outlineActor GetProperty] SetColor 0 0 0

# Add the actors to the renderer, set the background and size
#
ren1 AddActor outlineActor
ren1 AddActor isoActor
ren1 SetBackground 1 1 1
renWin SetSize 500 500
eval ren1 SetBackground $slate_grey

# render the image
#
iren SetUserMethod {wm deiconify .vtkInteract}
renWin Render

#renWin SetFileName "dicer.tcl.ppm"
#renWin SaveImageAsPPM

# prevent the tk window from showing up then start the event loop
wm withdraw .

