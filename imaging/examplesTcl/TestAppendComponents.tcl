# append multiple displaced spheres into an RGB image.
catch {load vtktcl}
source vtkImageInclude.tcl

# Image pipeline

vtkImageElipsoidSource sphere1
sphere1 SetCenter 95 100 0
sphere1 SetRadius 70 70 70

vtkImageElipsoidSource sphere2
sphere2 SetCenter 161 100 0
sphere2 SetRadius 70 70 70 

vtkImageElipsoidSource sphere3
sphere3 SetCenter 128 160 0
sphere3 SetRadius 70 70 70

vtkImageAppendComponents append1
append1 SetInput1 [sphere1 GetOutput]
append1 SetInput2 [sphere2 GetOutput]

vtkImageAppendComponents append2
append2 SetInput1 [sphere3 GetOutput]
append2 SetInput2 [append1 GetOutput]
append2 ReleaseDataFlagOff

vtkImageViewer viewer
viewer SetInput [append2 GetOutput]
viewer SetColorWindow 255
viewer SetColorLevel 128

# make interface
source WindowLevelInterface.tcl







