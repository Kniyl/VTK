package require vtk

# Image pipeline

vtkTIFFReader image1
  image1 SetFileName "$VTK_DATA_ROOT/Data/beach.tif"
  image1 Update

vtkStructuredPoints sp
eval sp SetDimensions [[image1 GetOutput] GetDimensions]
eval sp SetExtent [[image1 GetOutput] GetExtent]
sp SetScalarType [[image1 GetOutput] GetScalarType] 
sp SetNumberOfScalarComponents [[image1 GetOutput] GetNumberOfScalarComponents] 
[sp GetPointData] SetScalars [[[image1 GetOutput] GetPointData] GetScalars]

vtkImageLuminance luminance
  luminance SetInput sp

#
# write to the temp directory if possible, otherwise use .
#
set dir "."
if {[info commands rtTester] == "rtTester"}  {
   set dir [rtTester GetTempDirectory]
}

# make sure it is writeable first
if {[catch {set channel [open $dir/test.tmp w]}] == 0 } {
   close $channel
   file delete -force $dir/test.tmp

   vtkTIFFWriter tiff1
   tiff1 SetInput [image1 GetOutput]
   tiff1 SetFileName $dir/tiff1.tif
   
   vtkTIFFWriter tiff2
   tiff2 SetInput [luminance GetOutput]
   tiff2 SetFileName $dir/tiff2.tif
   
   vtkBMPWriter bmp1
   bmp1 SetInput [image1 GetOutput]
   bmp1 SetFileName $dir/bmp1.bmp
   
   vtkBMPWriter bmp2
   bmp2 SetInput [luminance GetOutput]
   bmp2 SetFileName $dir/bmp2.bmp
   
   vtkPNMWriter pnm1
   pnm1 SetInput [image1 GetOutput]
   pnm1 SetFileName $dir/pnm1.pnm
   
   vtkPNMWriter pnm2
   pnm2 SetInput [luminance GetOutput]
   pnm2 SetFileName $dir/pnm2.pnm

   vtkPostScriptWriter psw1
   psw1 SetInput [image1 GetOutput]
   psw1 SetFileName $dir/psw1.ps

   vtkPostScriptWriter psw2
   psw2 SetInput [luminance GetOutput]
   psw2 SetFileName $dir/psw2.ps

   vtkPNGWriter pngw1
   pngw1 SetInput [image1 GetOutput]
   pngw1 SetFileName $dir/pngw1.png

   vtkPNGWriter pngw2
   pngw2 SetInput [luminance GetOutput]
   pngw2 SetFileName $dir/pngw2.png

   vtkJPEGWriter jpgw1
   jpgw1 SetInput [image1 GetOutput]
   jpgw1 SetFileName $dir/jpgw1.jpg

   vtkJPEGWriter jpgw2
   jpgw2 SetInput [luminance GetOutput]
   jpgw2 SetFileName $dir/jpgw2.jpg

   tiff1 Write
   tiff2 Write
   bmp1 Write
   bmp2 Write
   pnm1 Write
   pnm2 Write
   psw1 Write
   psw2 Write
   pngw1 Write
   pngw2 Write
   jpgw1 Write
   jpgw2 Write

   file delete -force $dir/tiff1.tif
   file delete -force $dir/tiff2.tif
   file delete -force $dir/bmp1.bmp
   file delete -force $dir/bmp2.bmp
   file delete -force $dir/pnm1.pnm
   file delete -force $dir/pnm2.pnm
   file delete -force $dir/psw1.ps
   file delete -force $dir/psw2.ps
   file delete -force $dir/pngw1.png
   file delete -force $dir/pngw2.png
   file delete -force $dir/jpgw1.jpg
   file delete -force $dir/jpgw2.jpg
}

vtkImageViewer viewer
  viewer SetInput [luminance GetOutput]
  viewer SetColorWindow 255
  viewer SetColorLevel 127.5
  viewer Render

