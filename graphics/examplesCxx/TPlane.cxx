#include "vtk.h"

main ()
{
  vtkRenderer *ren1 = vtkRenderer::New();
  vtkRenderWindow *renWin = vtkRenderWindow::New();
    renWin->AddRenderer(ren1);
  vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
    iren->SetRenderWindow(renWin);

  // load the texture map
  vtkPNMReader *pnm = vtkPNMReader::New();
    pnm->SetFileName("../../../vtkdata/masonry.ppm");
  vtkTexture *atext = vtkTexture::New();
    atext->SetInput(pnm->GetOutput());
    atext->InterpolateOn();

  vtkPlaneSource *plane = vtkPlaneSource::New();
  vtkPolyDataMapper *planeMapper = vtkPolyDataMapper::New();
    planeMapper->SetInput(plane->GetOutput());
  vtkActor *planeActor = vtkActor::New();
    planeActor->SetMapper(planeMapper);
    planeActor->SetTexture(atext);

  ren1->AddActor(planeActor);
  ren1->SetBackground(0.2,0.3,0.4);
  renWin->SetSize(500,500);

  // interact with data
  renWin->Render();

  ren1->GetActiveCamera()->Zoom(1.4);
  renWin->Render();

  iren->Start();

  // Clean up
  ren1->Delete();
  renWin->Delete();
  iren->Delete();
  atext->Delete();
  pnm->Delete();
  plane->Delete();
  planeMapper->Delete();
  planeActor->Delete();
}
