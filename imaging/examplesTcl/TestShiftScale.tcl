catch {load vtktcl}
# Shift and scale an image (in that order)
# This filter is usefull for converting to a lower precision data type.


set sliceNumber 22

set VTK_FLOAT              1
set VTK_INT                2
set VTK_SHORT              3
set VTK_UNSIGNED_SHORT     4
set VTK_UNSIGNED_CHAR      5

set VTK_IMAGE_X_AXIS             0
set VTK_IMAGE_Y_AXIS             1
set VTK_IMAGE_Z_AXIS             2
set VTK_IMAGE_TIME_AXIS          3
set VTK_IMAGE_COMPONENT_AXIS     4


# Image pipeline

vtkImageReader reader
#reader DebugOn
[reader GetCache] ReleaseDataFlagOff
reader SetDataByteOrderToLittleEndian
reader SetDataExtent 0 255 0 255 1 93
reader SetFilePrefix "../../../vtkdata/fullHead/headsq"
reader SetDataMask 0x7fff

vtkImageShiftScale shiftScale
#shiftScale DebugOn
shiftScale SetInput [reader GetOutput]
shiftScale SetShift -3000.0
shiftScale SetScale -0.12

vtkImageViewer viewer
#viewer DebugOn
viewer SetInput [shiftScale GetOutput]
viewer SetZSlice $sliceNumber
viewer SetColorWindow 256
viewer SetColorLevel 128
viewer Render


#make interface
#

frame .slice
button .slice.up -text "Slice Up" -command SliceUp
button .slice.down -text "Slice Down" -command SliceDown

frame .wl
frame .wl.f1
label .wl.f1.windowLabel -text Window
scale .wl.f1.window -from 1 -to 512 -orient horizontal -command SetWindow \
  -variable window
frame .wl.f2
label .wl.f2.levelLabel -text Level
scale .wl.f2.level -from 1 -to 256 -orient horizontal -command SetLevel
checkbutton .wl.video -text "Inverse Video" -command SetInverseVideo


.wl.f1.window set 256
.wl.f2.level set 128


pack .slice .wl -side left
pack .slice.up .slice.down -side top
pack .wl.f1 .wl.f2 .wl.video -side top
pack .wl.f1.windowLabel .wl.f1.window -side left
pack .wl.f2.levelLabel .wl.f2.level -side left


proc SliceUp {} {
   global sliceNumber viewer
   if {$sliceNumber < 92} {set sliceNumber [expr $sliceNumber + 1]}
   puts $sliceNumber
   viewer SetZSlice $sliceNumber
   viewer Render
}

proc SliceDown {} {
   global sliceNumber viewer
   if {$sliceNumber > 0} {set sliceNumber [expr $sliceNumber - 1]}
   puts $sliceNumber
   viewer SetZSlice $sliceNumber
   viewer Render
}

proc SetWindow window {
   global viewer video
   if {$video} {
      viewer SetColorWindow [expr -$window]
   } else {
      viewer SetColorWindow $window
   }
   viewer Render
}

proc SetLevel level {
   global viewer
   viewer SetColorLevel $level
   viewer Render
}

proc SetInverseVideo {} {
   global viewer video window
   if {$video} {
      viewer SetColorWindow [expr -$window]
   } else {
      viewer SetColorWindow $window
   }
   viewer Render
}









