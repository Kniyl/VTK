catch {load vtktcl}
if { [catch {set VTK_TCL $env(VTK_TCL)}] != 0} { set VTK_TCL "../../examplesTcl" }
if { [catch {set VTK_DATA $env(VTK_DATA)}] != 0} { set VTK_DATA "../../../vtkdata" }


source vtkImageInclude.tcl
source TkImageViewerInteractor.tcl


vtkPNMReader reader
reader SetFileName "$VTK_DATA/vtks.pgm"

vtkImageFFT fft
fft SetFilteredAxes $VTK_IMAGE_X_AXIS $VTK_IMAGE_Y_AXIS
fft SetInput [reader GetOutput]

vtkImageMagnitude mag
mag SetInput [fft GetOutput]

vtkImageFourierCenter center
center SetInput [mag GetOutput]
center SetFilteredAxes $VTK_IMAGE_X_AXIS $VTK_IMAGE_Y_AXIS
 
vtkImageLogarithmicScale compress
compress SetInput [center GetOutput]
compress SetConstant 15

vtkImageViewer viewer1
viewer1 SetInput [reader GetOutput]
viewer1 SetColorWindow 160
viewer1 SetColorLevel 120


vtkImageViewer viewer
viewer SetInput [compress GetOutput]
viewer SetColorWindow 160
viewer SetColorLevel 120


# Create the GUI
wm withdraw .
toplevel .top 

frame .top.f1 

vtkTkImageViewerWidget .top.f1.r1 -width 320 -height 160 -iv viewer1
vtkTkImageViewerWidget .top.f1.r2 -width 320 -height 160 -iv viewer

button .top.btn  -text Quit -command exit

pack .top.f1.r1 .top.f1.r2 \
  -side top -padx 3 -pady 3 -expand t
pack .top.f1 -fill both -expand t
pack .top.btn -fill x

BindTkImageViewer .top.f1.r1 
BindTkImageViewer .top.f1.r2


