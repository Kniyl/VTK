#!/usr/bin/env python

# This example demonstrates the generation of a streamsurface.

import vtk
from vtk.util.misc import vtkGetDataRoot
VTK_DATA_ROOT = vtkGetDataRoot()

# Read the data and specify which scalars and vectors to read.
pl3d = vtk.vtkPLOT3DReader()
pl3d.SetXYZFileName(VTK_DATA_ROOT + "/Data/combxyz.bin")
pl3d.SetQFileName(VTK_DATA_ROOT + "/Data/combq.bin")
pl3d.SetScalarFunctionNumber(100)
pl3d.SetVectorFunctionNumber(202)
pl3d.Update()

# We use a rake to generate a series of streamline starting points
# scattered along a line. Each point will generate a streamline. These
# streamlines are then fed to the vtkRuledSurfaceFilter which stitches
# the lines together to form a surface.
rake = vtk.vtkLineSource()
rake.SetPoint1(15, -5, 32)
rake.SetPoint2(15, 5, 32)
rake.SetResolution(21)
rakeMapper = vtk.vtkPolyDataMapper()
rakeMapper.SetInput(rake.GetOutput())
rakeActor = vtk.vtkActor()
rakeActor.SetMapper(rakeMapper)

integ = vtk.vtkRungeKutta4()
sl = vtk.vtkStreamLine()
sl.SetInput(pl3d.GetOutput())
sl.SetSource(rake.GetOutput())
sl.SetIntegrator(integ)
sl.SetMaximumPropagationTime(0.1)
sl.SetIntegrationStepLength(0.1)
sl.SetIntegrationDirectionToBackward()
sl.SetStepLength(0.001)

# The ruled surface stiches together lines with triangle strips.
# Note the SetOnRatio method. It turns on every other strip that
# the filter generates (only when multiple lines are input).
scalarSurface = vtk.vtkRuledSurfaceFilter()
scalarSurface.SetInput(sl.GetOutput())
scalarSurface.SetOffset(0)
scalarSurface.SetOnRatio(2)
scalarSurface.PassLinesOn()
scalarSurface.SetRuledModeToPointWalk()
scalarSurface.SetDistanceFactor(30)
mapper = vtk.vtkPolyDataMapper()
mapper.SetInput(scalarSurface.GetOutput())
mapper.SetScalarRange(pl3d.GetOutput().GetScalarRange())
actor = vtk.vtkActor()
actor.SetMapper(mapper)

# Put an outline around for context.
outline = vtk.vtkStructuredGridOutlineFilter()
outline.SetInput(pl3d.GetOutput())
outlineMapper = vtk.vtkPolyDataMapper()
outlineMapper.SetInput(outline.GetOutput())
outlineActor = vtk.vtkActor()
outlineActor.SetMapper(outlineMapper)
outlineActor.GetProperty().SetColor(0, 0, 0)

# Now create the usual graphics stuff.
ren = vtk.vtkRenderer()
renWin = vtk.vtkRenderWindow()
renWin.AddRenderer(ren)
iren = vtk.vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin)

ren.AddActor(rakeActor)
ren.AddActor(actor)
ren.AddActor(outlineActor)
ren.SetBackground(1, 1, 1)

renWin.SetSize(300, 300)

iren.Initialize()
renWin.Render()
iren.Start()
