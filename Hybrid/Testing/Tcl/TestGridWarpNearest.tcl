package require vtktcl

# first, create an image to warp
vtkImageGridSource imageGrid
imageGrid SetGridSpacing 16 16 0
imageGrid SetGridOrigin 0 0 0
imageGrid SetDataExtent 0 255 0 255 0 0
imageGrid SetDataScalarTypeToUnsignedChar

vtkLookupTable table
table SetTableRange 0 1
table SetValueRange 1.0 0.0 
table SetSaturationRange 0.0 0.0 
table SetHueRange 0.0 0.0 
table SetAlphaRange 0.0 1.0
table Build

vtkImageMapToColors alpha
alpha SetInput [imageGrid GetOutput]
alpha SetLookupTable table

vtkBMPReader reader1
reader1 SetFileName "$VTK_DATA_ROOT/Data/masonry.bmp"

vtkImageBlend blend
blend SetInput 0 [reader1 GetOutput]
blend SetInput 1 [alpha GetOutput]

# next, create a ThinPlateSpline transform 

vtkPoints p1
p1 SetNumberOfPoints 8
p1 SetPoint 0 0 0 0
p1 SetPoint 1 0 255 0
p1 SetPoint 2 255 0 0
p1 SetPoint 3 255 255 0
p1 SetPoint 4 96 96 0
p1 SetPoint 5 96 159 0
p1 SetPoint 6 159 159 0
p1 SetPoint 7 159 96 0

vtkPoints p2
p2 SetNumberOfPoints 8
p2 SetPoint 0 0 0 0
p2 SetPoint 1 0 255 0
p2 SetPoint 2 255 0 0
p2 SetPoint 3 255 255 0
p2 SetPoint 4 96 159 0
p2 SetPoint 5 159 159 0
p2 SetPoint 6 159 96 0
p2 SetPoint 7 96 96 0

vtkThinPlateSplineTransform thinPlate
  thinPlate SetSourceLandmarks p2
  thinPlate SetTargetLandmarks p1
  thinPlate SetBasisToR2LogR

# convert the thin plate spline into a grid

# for nearest neighbor interpolation, the grid should precicely
# overlay the image you want to warp

vtkTransformToGrid transformToGrid
  transformToGrid SetInput thinPlate
  transformToGrid SetGridSpacing 1 1 1
  transformToGrid SetGridOrigin 0 0 0
  transformToGrid SetGridExtent 0 255 0 255 0 0

vtkGridTransform transform
  transform SetDisplacementGrid [transformToGrid GetOutput]
  transform SetInterpolationModeToNearestNeighbor
  # must lower the tolerance or it won't invert
  transform SetInverseTolerance 2.0
# you must invert the transform before passing it to vtkImageReslice
  transform Inverse

# apply the grid warp to the image

vtkImageReslice reslice
  reslice SetInput [blend GetOutput]
  reslice SetResliceTransform transform
  reslice SetInterpolationModeToLinear

# set the window/level to 255.0/127.5 to view full range
vtkImageViewer viewer
viewer SetInput [reslice GetOutput]
viewer SetColorWindow 255.0
viewer SetColorLevel 127.5
viewer SetZSlice 0
viewer Render



