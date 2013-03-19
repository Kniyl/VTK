import vtk, os, sys
from vtk.test import Testing

ss = vtk.vtkSphereSource() #make mesh to test with

af = vtk.vtkElevationFilter() #add some attributes
af.SetInputConnection(ss.GetOutputPort())

ef = vtk.vtkExtractEdges() #make lines to test
ef.SetInputConnection(af.GetOutputPort())

gf = vtk.vtkGlyph3D() #make verts to test
pts = vtk.vtkPoints()
pts.InsertNextPoint(0,0,0)
verts = vtk.vtkCellArray()
avert = vtk.vtkVertex()
avert.GetPointIds().SetId(0, 0)
verts.InsertNextCell(avert)
onevertglyph = vtk.vtkPolyData()
onevertglyph.SetPoints(pts)
onevertglyph.SetVerts(verts)
gf.SetSourceData(onevertglyph)
gf.SetInputConnection(af.GetOutputPort())

testwrites = ["points","lines","mesh"]
failed = False
for datasetString in testwrites:
  if datasetString == "points":
    toshow=gf
  elif datasetString == "lines":
    toshow = ef
  else:
    toshow = af
  gw = vtk.vtkGeoJSONWriter()
  gw.SetInputConnection(toshow.GetOutputPort())
  gw.SetFileName("sphere.gjson")
  gw.Write()
  if os.path.exists("sphere.gjson") and os.path.isfile("sphere.gjson"):
    os.remove("sphere.gjson")
  else:
    print "Failed to write " + datasetString + " to file"
    failed = True
  gw.WriteToOutputStringOn()
  gw.Write()
  gj = "['"+str(gw.RegisterAndGetOutputString()).replace('\n','')+"']"
  if len(gj) <= 1000:
    print "Failed to write " + datasetString + " to buffer"
    failed = True

if failed:
  sys.exit(1)
sys.exit(0)
