catch {load vtktcl}
# get the interactor ui
source ../../examplesTcl/vtkInt.tcl

# Simple volume rendering example.
vtkSLCReader reader
    reader SetFileName "../../../vtkdata/poship.slc"

# Create transfer functions for opacity and color
vtkPiecewiseFunction opacityTransferFunction
    opacityTransferFunction AddPoint  20   0.0
    opacityTransferFunction AddPoint  255  0.2

vtkColorTransferFunction colorTransferFunction
    colorTransferFunction AddRGBPoint      0.0 0.0 0.0 0.0
    colorTransferFunction AddRGBPoint     64.0 1.0 0.0 0.0
    colorTransferFunction AddRGBPoint    128.0 0.0 0.0 1.0
    colorTransferFunction AddRGBPoint    192.0 0.0 1.0 0.0
    colorTransferFunction AddRGBPoint    255.0 0.0 0.2 0.0

# Create property, mapper, volume 
vtkVolumeProperty volumeProperty
    volumeProperty SetColor colorTransferFunction
    volumeProperty SetScalarOpacity opacityTransferFunction

vtkVolumeTextureMapper2D volumeMapper
    volumeMapper SetInput [reader GetOutput]

vtkVolume volume
    volume SetMapper volumeMapper
    volume SetProperty volumeProperty

# Okay now the graphics stuff
vtkRenderer ren1
vtkRenderWindow renWin
    renWin AddRenderer ren1
    renWin SetSize 256 256
vtkRenderWindowInteractor iren
    iren SetRenderWindow renWin

ren1 AddVolume volume
ren1 SetBackground 0.1 0.2 0.4
renWin Render

proc TkCheckAbort {} {
  set foo [renWin GetEventPending]
  if {$foo != 0} {renWin SetAbortRender 1}
}
renWin SetAbortCheckMethod {TkCheckAbort}

iren SetUserMethod {wm deiconify .vtkInteract}
iren Initialize

#renWin SetFileName "volTexSimple.tcl.ppm"
#renWin SaveImageAsPPM

wm withdraw .











