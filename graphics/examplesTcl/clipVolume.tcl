# Example demonstrates how to generate a 3D tetrahedra mesh from a volume
#
catch {load vtktcl}
if { [catch {set VTK_TCL $env(VTK_TCL)}] != 0} { set VTK_TCL "../../examplesTcl" }
if { [catch {set VTK_DATA $env(VTK_DATA)}] != 0} { set VTK_DATA "../../../vtkdata" }

source $VTK_TCL/vtkInt.tcl

# Quadric definition
vtkQuadric quadric
  quadric SetCoefficients .5 1 .2 0 .1 0 0 .2 0 0

vtkSampleFunction sample
  sample SetSampleDimensions 20 20 20
  sample SetImplicitFunction quadric
  sample ComputeNormalsOff
    
# Generate tetrahedral mesh
vtkClipVolume clip
  clip SetInput [sample GetOutput]
  clip SetValue 1.0
  clip GenerateClippedOutputOff

vtkDataSetMapper clipMapper
  clipMapper SetInput [clip GetOutput]
  clipMapper ScalarVisibilityOff

vtkActor clipActor
  clipActor SetMapper clipMapper
  [clipActor GetProperty] SetColor .8 .4 .4

# Create outline
vtkOutlineFilter outline
  outline SetInput [clip GetInput]

vtkPolyDataMapper outlineMapper
  outlineMapper SetInput [outline GetOutput]

vtkActor outlineActor
  outlineActor SetMapper outlineMapper
  eval [outlineActor GetProperty] SetColor 0 0 0

# Define graphics objects
vtkRenderer ren1
vtkRenderWindow renWin
    renWin AddRenderer ren1
vtkRenderWindowInteractor iren
    iren SetRenderWindow renWin

ren1 SetBackground 1 1 1
ren1 AddActor clipActor
ren1 AddActor outlineActor

iren SetUserMethod {wm deiconify .vtkInteract}
iren Initialize

#renWin SetFileName clipVolume.tcl.ppm
#renWin SaveImageAsPPM

wm withdraw .
