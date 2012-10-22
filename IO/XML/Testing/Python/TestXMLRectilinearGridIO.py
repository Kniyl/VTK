#!/usr/bin/env python

file0 = rgFile0.vtr
file1 = rgFile1.vtr
file2 = rgFile2.vtr
# read in some grid data
gridReader = vtk.vtkRectilinearGridReader()
gridReader.SetFileName("" + str(VTK_DATA_ROOT) + "/Data/RectGrid2.vtk")
gridReader.Update()
# extract to reduce extents of grid
extract = vtk.vtkExtractRectilinearGrid()
extract.SetInputConnection(gridReader.GetOutputPort())
extract.SetVOI(0,23,0,32,0,10)
extract.Update()
# write just a piece (extracted piece) as well as the whole thing
rgWriter = vtk.vtkXMLRectilinearGridWriter()
rgWriter.SetFileName(file0)
rgWriter.SetInputConnection(extract.GetOutputPort())
rgWriter.SetDataModeToAscii()
rgWriter.Write()
rgWriter.SetFileName(file1)
rgWriter.SetInputConnection(gridReader.GetOutputPort())
rgWriter.SetDataModeToAppended()
rgWriter.SetNumberOfPieces(2)
rgWriter.Write()
rgWriter.SetFileName(file2)
rgWriter.SetDataModeToBinary()
rgWriter.SetWriteExtent(3,46,6,32,1,5)
rgWriter.SetCompressor("")
if (rgWriter.GetByteOrder()):
    rgWriter.SetByteOrder(0)
    pass
else:
    rgWriter.SetByteOrder(1)
    pass
rgWriter.Write()
# read the extracted grid
reader = vtk.vtkXMLRectilinearGridReader()
reader.SetFileName(file0)
reader.WholeSlicesOff()
reader.Update()
rg0 = vtk.vtkRectilinearGrid()
rg0.DeepCopy(reader.GetOutput())
mapper0 = vtk.vtkDataSetMapper()
mapper0.SetInputData(rg0)
actor0 = vtk.vtkActor()
actor0.SetMapper(mapper0)
# read the whole grid
reader.SetFileName(file1)
reader.WholeSlicesOn()
reader.Update()
rg1 = vtk.vtkRectilinearGrid()
rg1.DeepCopy(reader.GetOutput())
mapper1 = vtk.vtkDataSetMapper()
mapper1.SetInputData(rg1)
actor1 = vtk.vtkActor()
actor1.SetMapper(mapper1)
actor1.SetPosition(-1.5,3,0)
# read the partially written grid
reader.SetFileName(file2)
reader.Update()
mapper2 = vtk.vtkDataSetMapper()
mapper2.SetInputConnection(reader.GetOutputPort())
actor2 = vtk.vtkActor()
actor2.SetMapper(mapper2)
actor2.SetPosition(1.5,3,0)
# Create the RenderWindow, Renderer and both Actors
#
ren1 = vtk.vtkRenderer()
renWin = vtk.vtkRenderWindow()
renWin.AddRenderer(ren1)
iren = vtk.vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin)
# Add the actors to the renderer, set the background and size
#
ren1.AddActor(actor0)
ren1.AddActor(actor1)
ren1.AddActor(actor2)
renWin.SetSize(300,300)
iren.Initialize()
renWin.Render()
# prevent the tk window from showing up then start the event loop
file.delete("-force", file0)
file.delete("-force", file1)
file.delete("-force", file2)
# --- end of script --
