#!/usr/bin/env python
import vtk
from vtk.test import Testing
from vtk.util.misc import vtkGetDataRoot
VTK_DATA_ROOT = vtkGetDataRoot()

# first, create an image to warp
imageGrid = vtk.vtkImageGridSource()
imageGrid.SetGridSpacing(16,16,0)
imageGrid.SetGridOrigin(0,0,0)
imageGrid.SetDataExtent(0,255,0,255,0,0)
imageGrid.SetDataScalarTypeToUnsignedChar()
table = vtk.vtkLookupTable()
table.SetTableRange(0,1)
table.SetValueRange(1.0,0.0)
table.SetSaturationRange(0.0,0.0)
table.SetHueRange(0.0,0.0)
table.SetAlphaRange(0.0,1.0)
table.Build()
alpha = vtk.vtkImageMapToColors()
alpha.SetInputConnection(imageGrid.GetOutputPort())
alpha.SetLookupTable(table)
reader1 = vtk.vtkBMPReader()
reader1.SetFileName("" + str(VTK_DATA_ROOT) + "/Data/masonry.bmp")
blend = vtk.vtkImageBlend()
blend.AddInputConnection(reader1.GetOutputPort())
blend.AddInputConnection(alpha.GetOutputPort())
# next, create a ThinPlateSpline transform
p1 = vtk.vtkPoints()
p1.SetNumberOfPoints(8)
p1.SetPoint(0,0,0,0)
p1.SetPoint(1,0,255,0)
p1.SetPoint(2,255,0,0)
p1.SetPoint(3,255,255,0)
p1.SetPoint(4,96,96,0)
p1.SetPoint(5,96,159,0)
p1.SetPoint(6,159,159,0)
p1.SetPoint(7,159,96,0)
p2 = vtk.vtkPoints()
p2.SetNumberOfPoints(8)
p2.SetPoint(0,0,0,0)
p2.SetPoint(1,0,255,0)
p2.SetPoint(2,255,0,0)
p2.SetPoint(3,255,255,0)
p2.SetPoint(4,96,159,0)
p2.SetPoint(5,159,159,0)
p2.SetPoint(6,159,96,0)
p2.SetPoint(7,96,96,0)
thinPlate = vtk.vtkThinPlateSplineTransform()
thinPlate.SetSourceLandmarks(p2)
thinPlate.SetTargetLandmarks(p1)
thinPlate.SetBasisToR2LogR()
# convert the thin plate spline into a grid
# for nearest neighbor interpolation, the grid should precicely
# overlay the image you want to warp
transformToGrid = vtk.vtkTransformToGrid()
transformToGrid.SetInput(thinPlate)
transformToGrid.SetGridSpacing(1,1,1)
transformToGrid.SetGridOrigin(0,0,0)
transformToGrid.SetGridExtent(0,255,0,255,0,0)
transformToGrid.Update()
transform = vtk.vtkGridTransform()
transform.SetDisplacementGridConnection(transformToGrid.GetOutputPort())
transform.SetInterpolationModeToNearestNeighbor()
# must lower the tolerance or it won't invert
transform.SetInverseTolerance(2.0)
# you must invert the transform before passing it to vtkImageReslice
transform.Inverse()
# apply the grid warp to the image
reslice = vtk.vtkImageReslice()
reslice.SetInputConnection(blend.GetOutputPort())
reslice.SetResliceTransform(transform)
reslice.SetInterpolationModeToLinear()
# set the window/level to 255.0/127.5 to view full range
viewer = vtk.vtkImageViewer()
viewer.SetInputConnection(reslice.GetOutputPort())
viewer.SetColorWindow(255.0)
viewer.SetColorLevel(127.5)
viewer.SetZSlice(0)
viewer.Render()
# --- end of script --
