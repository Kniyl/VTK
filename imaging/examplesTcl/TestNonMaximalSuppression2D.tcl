# This script is for testing the 2dNonMaximalSuppressionFilter.
# The filter works exclusively on the output of the gradient filter,
# because it expects a phase (in radians) in the second component.
# The effect is to pick the peaks of the gradient creating thin lines.


#source ../../examples/tcl/colors.tcl
#source ../../examples/tcl/vtkInt.tcl

set sliceNumber 45

set VTK_IMAGE_FLOAT              1
set VTK_IMAGE_INT                2
set VTK_IMAGE_SHORT              3
set VTK_IMAGE_UNSIGNED_SHORT     4
set VTK_IMAGE_UNSIGNED_CHAR      5

set VTK_IMAGE_X_AXIS             0
set VTK_IMAGE_Y_AXIS             1
set VTK_IMAGE_Z_AXIS             2
set VTK_IMAGE_TIME_AXIS          3

# First create the render master
vtkRenderMaster rm;





# Image pipeline
vtkImage4dShortReader reader;
#reader DebugOn
[reader GetCache] ReleaseDataFlagOff;
reader SwapBytesOn;
reader SetSize 256 256 94 1;
reader SetFileRoot "../../data/fullHead/headsq.%d"
reader SetPixelMask 0x7fff;

vtkImage2dGradientFilter gradient;
gradient SetInput [reader GetOutput];

vtkImage2dNonMaximalSuppressionFilter suppress;
suppress SetInput [gradient GetOutput];
suppress ReleaseDataFlagOff;

vtkImageXViewer viewer;
#viewer DebugOn;
viewer SetInput [suppress GetOutput];
viewer SetDefaultCoordinate2 $sliceNumber;
viewer SetColorWindow 100
viewer SetColorLevel 40
viewer Render;


#make interface
#
frame .slice
button .slice.up -text "Slice Up" -command SliceUp
button .slice.down -text "Slice Down" -command SliceDown

frame .wl
frame .wl.f1;
label .wl.f1.windowLabel -text Window;
scale .wl.f1.window -from 1 -to 500 -orient horizontal -command SetWindow
frame .wl.f2;
label .wl.f2.levelLabel -text Level;
scale .wl.f2.level -from -200 -to 200 -orient horizontal -command SetLevel
checkbutton .wl.video -text "Inverse Video" -variable inverseVideo -command SetInverseVideo


.wl.f1.window set 100
.wl.f2.level set 40


pack .slice .wl -side left
pack .slice.up .slice.down -side top
pack .wl.f1 .wl.f2 .wl.video -side top
pack .wl.f1.windowLabel .wl.f1.window -side left
pack .wl.f2.levelLabel .wl.f2.level -side left

#$renWin SetTkWindow .renwin

proc SliceUp {} {
        global sliceNumber viewer
        if {$sliceNumber < 93} {set sliceNumber [expr $sliceNumber + 1]}
        viewer SetDefaultCoordinate2 $sliceNumber;
	viewer Render
        puts $sliceNumber;
}

proc SliceDown {} {
        global sliceNumber tform renWin
        if {$sliceNumber > 0} {set sliceNumber [expr $sliceNumber - 1]}
        viewer SetDefaultCoordinate2 $sliceNumber;
	viewer Render
        puts $sliceNumber;
}

proc SetWindow window {
	global viewer
        viewer SetColorWindow $window;
        viewer Render;  
}
proc SetLevel level {
	global viewer
        viewer SetColorLevel $level;
        viewer Render;  
}
proc SetInverseVideo {} {
}







