catch {load vtktcl}
source ../../examplesTcl/vtkInt.tcl


vtkSLCReader reader
    reader SetFileName ../../../vtkdata/poship.slc

vtkSLCReader reader2
    reader2 SetFileName ../../../vtkdata/neghip.slc

vtkPiecewiseFunction opacityTransferFunction
    opacityTransferFunction AddPoint  20   0.0
    opacityTransferFunction AddPoint  255  0.3

vtkColorTransferFunction colorTransferFunction
    colorTransferFunction AddRedPoint      0.0 0.0
    colorTransferFunction AddRedPoint     64.0 1.0
    colorTransferFunction AddRedPoint    128.0 0.0
    colorTransferFunction AddRedPoint    255.0 0.0
    colorTransferFunction AddBluePoint    0.0 0.0
    colorTransferFunction AddBluePoint   64.0 0.0
    colorTransferFunction AddBluePoint  128.0 1.0
    colorTransferFunction AddBluePoint  192.0 0.0
    colorTransferFunction AddBluePoint  255.0 0.0
    colorTransferFunction AddGreenPoint     0.0 0.0
    colorTransferFunction AddGreenPoint   128.0 0.0
    colorTransferFunction AddGreenPoint   192.0 1.0
    colorTransferFunction AddGreenPoint   255.0 0.2

vtkVolumeProperty volumeProperty
    volumeProperty SetColor colorTransferFunction
    volumeProperty SetScalarOpacity opacityTransferFunction
    volumeProperty SetInterpolationTypeToLinear
    volumeProperty ShadeOn

vtkVolumeRayCastMIPFunction  MIPFunction

vtkVolumeRayCastMapper volumeMapper
    volumeMapper SetScalarInput [reader GetOutput]
    volumeMapper SetVolumeRayCastFunction MIPFunction
    volumeMapper SetSampleDistance 0.25

vtkVolume volume
    volume SetVolumeMapper volumeMapper
    volume SetVolumeProperty volumeProperty

vtkContourFilter contour
  contour SetInput [reader2 GetOutput]
  contour SetValue 0 128.0

vtkPolyDataMapper neghip_mapper
  neghip_mapper SetInput [contour GetOutput]
  neghip_mapper ScalarVisibilityOff

vtkActor neghip
  neghip SetMapper neghip_mapper
  [neghip GetProperty] SetColor 0.8 0.2 0.8
  [neghip GetProperty] SetAmbient 0.1
  [neghip GetProperty] SetDiffuse 0.6
  [neghip GetProperty] SetSpecular 0.4

# Okay now the graphics stuff
vtkRenderer ren1
    ren1 SetBackground 0.1 0.2 0.4
vtkRenderWindow renWin
    renWin AddRenderer ren1
    renWin SetSize 256 256
vtkRenderWindowInteractor iren
    iren SetRenderWindow renWin

ren1 AddActor neghip
ren1 SetBackground 0 0 0 

[ren1 GetActiveCamera] SetPosition 162.764 30.8946 116.029
[ren1 GetActiveCamera] SetFocalPoint 32.868 31.5566 31.9246
[ren1 GetActiveCamera] SetViewUp -0.00727828 0.999791 0.0191114
[ren1 GetActiveCamera] SetViewPlaneNormal 0.839404 -0.00427837 0.543492
[ren1 GetActiveCamera] SetClippingRange 15.4748 773.74

ren1 AddVolume volume
renWin SetSize 200 200
renWin Render

iren SetUserMethod {wm deiconify .vtkInteract}
iren SetDesiredUpdateRate 1
iren Initialize

#renWin SetFileName "VolMIP.tcl.ppm"
#renWin SaveImageAsPPM

wm withdraw .



