#
#   Clip Actor with Spherical Lens
#

# get the interactor ui
source ../../examplesTcl/vtkInt.tcl
source ../../examplesTcl/colors.tcl

proc mkname {a b} {return $a$b}

# proc to make actors
# create pipeline
proc MakeActor { name r g b} {
#
    set filename  [eval mkname $name .vtk]
    set reader  [eval mkname $name PolyDataReader]
    vtkPolyDataReader $reader
       $reader SetFileName $filename
    set mapper [eval mkname $name PolyDataMapper]
    vtkPolyDataMapper $mapper
        $mapper SetInput [$reader GetOutput]
        $mapper ScalarVisibilityOff
        $mapper ImmediateModeRenderingOn
    set actor [ eval mkname $name Actor]
    vtkActor $actor
        $actor SetMapper $mapper
        eval [$actor GetProperty] SetDiffuseColor $r $g $b
        eval [$actor GetProperty] SetSpecularPower 50
        eval [$actor GetProperty] SetSpecular .5
        eval [$actor GetProperty] SetDiffuse .8
    return $actor
}

# Now create the RenderWindow, Renderer and Interactor
#
vtkRenderer ren1
vtkRenderWindow renWin
  renWin AddRenderer ren1
vtkRenderWindowInteractor iren
  iren SetRenderWindow renWin

# Add the actors to the renderer, set the background and size
#

ren1 AddActor [eval MakeActor ../../../vtkdata/skin $flesh]

set x 0
set y 70
set z -70
set r 30

#
# make a base plane
#
vtkPlaneSource basePlane
  basePlane SetOrigin -102.4 139.5 -102.4
  basePlane SetPoint1 101.6 139.5 -102.4
  basePlane SetPoint2 -102.4 139.5 101.6
vtkPolyDataMapper baseMapper
  baseMapper SetInput [basePlane GetOutput]
vtkActor base
  base SetMapper baseMapper
ren1 AddActor base

#
# make the geometry for a lens
#
vtkSphereSource lensSource
  lensSource SetRadius $r
  lensSource SetCenter $x $y $z
  lensSource SetThetaResolution 256
  lensSource SetPhiResolution 256

vtkPolyDataMapper lensGeometryMapper
  lensGeometryMapper SetInput [lensSource GetOutput]
  lensGeometryMapper ImmediateModeRenderingOn

vtkActor lensGeometry
  lensGeometry SetMapper lensGeometryMapper
  lensGeometry VisibilityOff

#
# read the volume
#
set RESOLUTION 256
set START_SLICE 1
set END_SLICE 93
set PIXEL_SIZE .8
set origin [expr ( $RESOLUTION / 2.0 ) * $PIXEL_SIZE * -1.0]
set SLICE_ORDER si
source ../../examplesTcl/frog/SliceOrder.tcl

vtkVolume16Reader reader
  eval reader SetDataDimensions $RESOLUTION $RESOLUTION
  eval reader SetFilePrefix ../../../vtkdata/fullHead/headsq
  eval reader SetDataSpacing $PIXEL_SIZE $PIXEL_SIZE 1.5
  eval reader SetDataOrigin $origin $origin 1.5
  eval reader SetImageRange $START_SLICE $END_SLICE
  reader SetTransform si
  reader SetHeaderSize 0
  reader SetDataMask 0x7fff
  reader SetDataByteOrderToLittleEndian
  reader Update

vtkImplicitVolume aVolumeModel
  aVolumeModel SetVolume [reader GetOutput]
  aVolumeModel SetOutValue 0

#
# clip the lens geometry
#
vtkClipPolyData lensClipper
  lensClipper SetInput [lensSource GetOutput]
  lensClipper SetClipFunction aVolumeModel
  lensClipper SetValue 600.5
  lensClipper GenerateClipScalarsOn
  lensClipper GenerateClippedOutputOff
  lensClipper Update

vtkWindowLevelLookupTable wlLut 
  wlLut SetWindow 1000
  wlLut SetLevel 1200
  wlLut SetTableRange 0 2047
  wlLut Build

vtkPolyDataMapper lensMapper
  lensMapper SetInput [lensClipper GetOutput]
  lensMapper SetScalarRange 0 2047
  lensMapper SetLookupTable wlLut
  lensMapper ScalarVisibilityOn

vtkActor lens
  lens SetMapper lensMapper

#
# clip the surface geometry with the lens function
#
vtkSphere lensFunction
  lensFunction SetCenter $x $y $z
  lensFunction SetRadius $r

vtkClipPolyData surfaceClipper
  surfaceClipper SetInput [../../../vtkdata/skinPolyDataReader GetOutput]
  surfaceClipper SetClipFunction lensFunction
  surfaceClipper GenerateClippedOutputOn
  surfaceClipper GenerateClipScalarsOn
  surfaceClipper InsideOutOn
  surfaceClipper Update

vtkPolyDataMapper surfaceMapper
  surfaceMapper SetInput [surfaceClipper GetOutput]
  surfaceMapper SetScalarRange -100 100
  surfaceMapper ScalarVisibilityOff
  surfaceMapper ImmediateModeRenderingOn

vtkActor clippedSurface
  clippedSurface SetMapper surfaceMapper
eval [clippedSurface GetProperty] SetDiffuseColor $banana
  [clippedSurface GetProperty] SetSpecular .4
  [clippedSurface GetProperty] SetSpecularPower 30
  [clippedSurface GetProperty] SetOpacity .5
     
vtkPolyDataMapper insideSurfaceMapper
  insideSurfaceMapper SetInput [surfaceClipper GetClippedOutput]
  insideSurfaceMapper ScalarVisibilityOff
  insideSurfaceMapper ImmediateModeRenderingOn

../../../vtkdata/skinActor SetMapper insideSurfaceMapper
../../../vtkdata/skinActor VisibilityOn

#
# set up volume rendering
#
vtkPiecewiseFunction tfun
  tfun AddPoint   70.0  0.0
  tfun AddPoint  599.0  0
  tfun AddPoint  600.0 0
  tfun AddPoint  1195.0 0
  tfun AddPoint  1200 .2
  tfun AddPoint  1300 .3
  tfun AddPoint  2000 .3
  tfun AddPoint  4095.0  1.0

vtkColorTransferFunction ctfun
  ctfun AddRedPoint     0.0 0.5
  ctfun AddRedPoint    600.0 1.0
  ctfun AddRedPoint   4095.0 0.5

  ctfun AddGreenPoint  600.0 0.5
  ctfun AddGreenPoint 1280.0 .2
  ctfun AddGreenPoint 4095.0 0.5

  ctfun AddBluePoint  600.0 0.5
  ctfun AddBluePoint  1960.0 .1
  ctfun AddBluePoint  4095.0 0.5

vtkVolumeRayCastCompositeFunction  compositeFunction

vtkProjectedPolyDataRayBounder raybounder
  raybounder SetPolyData [ lensSource GetOutput]

vtkVolumeRayCastMapper volumeMapper
  volumeMapper SetInput [reader GetOutput]
  volumeMapper SetVolumeRayCastFunction compositeFunction
  volumeMapper SetRayBounder raybounder

vtkVolumeProperty volumeProperty
  volumeProperty SetColor ctfun
  volumeProperty SetScalarOpacity tfun
  volumeProperty SetInterpolationTypeToLinear
  volumeProperty ShadeOn

vtkVolume newvol
  newvol SetMapper volumeMapper
  newvol SetProperty volumeProperty

ren1 AddVolume newvol

ren1 AddActor lens
lens PickableOff
lens VisibilityOn

#ren1 AddActor clippedSurface
clippedSurface VisibilityOff
ren1 AddActor lensGeometry

lensGeometry VisibilityOff

ren1 SetBackground 0.2 0.3 0.4

renWin SetSize 320 240
[ren1 GetActiveCamera] SetViewUp 0 -1 0
[ren1 GetActiveCamera] Azimuth 230
[ren1 GetActiveCamera] Elevation 30
[ren1 GetActiveCamera] Dolly 1.75
ren1 ResetCameraClippingRange

iren Initialize

#
# Clip with spherical lens
#
proc Clip {x y z r} {
  lensSource SetCenter $x $y $z
  lensSource SetRadius $r
  lensFunction SetCenter $x $y $z
  lensFunction SetRadius $r
puts "lensClipper [expr [lindex [time {lensClipper Update;} 1] 0] / 1000000.0] seconds"
puts "surfaceClipper [expr [lindex [time {surfaceClipper Update;} 1] 0] / 1000000.0] seconds"
  renWin Render
}

proc PickAndClip {} {
  global r
  eval Clip [[iren GetPicker] GetPickPosition] $r
}

iren SetDesiredUpdateRate .5
iren SetEndPickMethod  {PickAndClip}
proc PickAndPrint {} {
  puts "eval [[iren GetPicker] GetPickPosition]"
}

# render the image
#
iren SetUserMethod {wm deiconify .vtkInteract}
#renWin SetFileName surfVol.tcl.ppm
#renWin SaveImageAsPPM
# prevent the tk window from showing up then start the event loop
wm withdraw .
