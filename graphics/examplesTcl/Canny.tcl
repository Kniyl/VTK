catch {load vtktcl}
if { [catch {set VTK_TCL $env(VTK_TCL)}] != 0} { set VTK_TCL "../../examplesTcl" }
if { [catch {set VTK_DATA $env(VTK_DATA)}] != 0} { set VTK_DATA "../../../vtkdata" }


# get the interactor ui
source $VTK_TCL/vtkInt.tcl

# Create the RenderWindow, Renderer and both Actors
vtkRenderer ren1
vtkRenderWindow renWin
    renWin AddRenderer ren1
vtkRenderWindowInteractor iren
    iren SetRenderWindow renWin

# load in the texture map
#
vtkBMPReader imageIn
imageIn SetFileName "$VTK_DATA/beach.bmp"

vtkImageLuminance il
il SetInput [imageIn GetOutput]

vtkImageCast ic
ic SetOutputScalarTypeToFloat
ic SetInput [il GetOutput]

# smooth the image
vtkImageGaussianSmooth gs
gs SetInput [ic GetOutput]
gs SetDimensionality 2
gs SetRadiusFactors 1 1 0

# gradient the image
vtkImageGradient imgGradient;
imgGradient SetInput [gs GetOutput];
imgGradient SetDimensionality 2

vtkImageMagnitude imgMagnitude
imgMagnitude SetInput [imgGradient GetOutput]

# non maximum suppression
vtkImageNonMaximumSuppression nonMax;
nonMax SetVectorInput [imgGradient GetOutput];
nonMax SetMagnitudeInput [imgMagnitude GetOutput];
nonMax SetDimensionality 2

vtkImageConstantPad pad
pad SetInput [imgGradient GetOutput]
pad SetOutputNumberOfScalarComponents 3
pad SetConstant 0

vtkImageToStructuredPoints i2sp1
i2sp1 SetInput [nonMax GetOutput]
i2sp1 SetVectorInput [pad GetOutput]

# link edgles
vtkLinkEdgels imgLink;
imgLink SetInput [i2sp1 GetOutput];
imgLink SetGradientThreshold 2;

# threshold links
vtkThreshold thresholdEdgels;
thresholdEdgels SetInput [imgLink GetOutput];
thresholdEdgels ThresholdByUpper 10;
thresholdEdgels AllScalarsOff

vtkGeometryFilter gf
gf SetInput [thresholdEdgels GetOutput]

vtkImageToStructuredPoints i2sp
i2sp SetInput [imgMagnitude GetOutput]
i2sp SetVectorInput [pad GetOutput]

# subpixel them
vtkSubPixelPositionEdgels spe;
spe SetInput [gf GetOutput];
spe SetGradMaps [i2sp GetOutput];

vtkStripper strip
strip SetInput [spe GetOutput]

vtkPolyDataMapper dsm
dsm SetInput [strip GetOutput]
dsm SetScalarRange 0 70

vtkActor planeActor
planeActor SetMapper dsm
[planeActor GetProperty] SetAmbient .5
[planeActor GetProperty] SetDiffuse 1.0

# Add the actors to the renderer, set the background and size
ren1 AddActor planeActor
ren1 SetBackground 0.1 0.2 0.4
renWin SetSize 500 500

# render the image
iren Initialize
iren SetUserMethod {wm deiconify .vtkInteract}
renWin Render
[ren1 GetActiveCamera] Zoom 1.4
renWin Render

# prevent the tk window from showing up then start the event loop
wm withdraw .





