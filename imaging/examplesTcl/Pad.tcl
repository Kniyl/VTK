catch {load vtktcl}
if { [catch {set VTK_TCL $env(VTK_TCL)}] != 0} { set VTK_TCL "../../examplesTcl" }
if { [catch {set VTK_DATA $env(VTK_DATA)}] != 0} { set VTK_DATA "../../../vtkdata" }

# Make an image larger by repeating the data.  Tile.


source vtkImageInclude.tcl
source TkImageViewerInteractor.tcl



# Image pipeline

vtkImageReader reader
reader SetDataByteOrderToLittleEndian
reader SetDataExtent 0 255 0 255 1 94
reader SetFilePrefix "$VTK_DATA/fullHead/headsq"
reader SetDataMask 0x7fff
#reader DebugOn

vtkImageMirrorPad pad1
pad1 SetInput [reader GetOutput]
pad1 SetOutputWholeExtent -127 383 -127 383 0 93

vtkImageConstantPad pad2
pad2 SetInput [reader GetOutput]
pad2 SetOutputWholeExtent -127 383 -127 383 0 93
pad2 SetConstant 800

vtkImageViewer viewer
viewer SetInput [pad1 GetOutput]
viewer SetZSlice 22
viewer SetColorWindow 2000
viewer SetColorLevel 1000
[viewer GetActor2D] SetDisplayPosition 127 127

vtkImageViewer viewer2
viewer2 SetInput [pad2 GetOutput]
viewer2 SetZSlice 22
viewer2 SetColorWindow 2000
viewer2 SetColorLevel 1000
[viewer2 GetActor2D] SetDisplayPosition 127 127


# Create the GUI
wm withdraw .
toplevel .top 

frame .top.f1 

vtkTkImageViewerWidget .top.f1.r1 -width 360 -height 512 -iv viewer
vtkTkImageViewerWidget .top.f1.r2 -width 360 -height 512 -iv viewer2

button .top.btn  -text Quit -command exit

pack .top.f1.r1 .top.f1.r2 \
  -side left -padx 3 -pady 3 -expand t
pack .top.f1  -fill both -expand t
pack .top.btn -fill x


BindTkImageViewer .top.f1.r1 
BindTkImageViewer .top.f1.r2



