// This program test the ports by setting up a simple pipeline.

#include "mpi.h"
#include "vtkImageGaussianSource.h"
#include "vtkImageEllipsoidSource.h"
#include "vtkImageToStructuredPoints.h"
#include "vtkUpStreamPort.h"
#include "vtkDownStreamPort.h"
#include "vtkTexture.h"
#include "vtkPlaneSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"


VTK_THREAD_RETURN_TYPE process_a( void *vtkNotUsed(arg) )
{
  vtkMultiProcessController *controller;
  vtkImageGaussianSource *source = vtkImageGaussianSource::New();
  vtkImageEllipsoidSource *ellipse = vtkImageEllipsoidSource::New();
  vtkUpStreamPort *upStreamPort = vtkUpStreamPort::New();
  int myid;
  
  controller = vtkMultiProcessController::RegisterAndGetGlobalController(NULL);
  myid = controller->GetLocalProcessId();
  
  // Set up the pipeline source.
  source->SetCenter(128.0, 128.0, 0.0);
  source->SetMaximum(2.0);
  source->SetStandardDeviation(50.0);

  ellipse->SetCenter(128.0, 128.0, 0.0);
  ellipse->SetRadius(50.0, 70.0, 1.0);
  
  vtkImageToStructuredPoints *sp = vtkImageToStructuredPoints::New();
  sp->SetInput(source->GetOutput());

  upStreamPort->SetInput((vtkImageData*)(sp->GetOutput()));
  upStreamPort->SetTag(999);
  
  // wait for the call back to execute.
  upStreamPort->WaitForUpdate();
  
  source->Delete();
  upStreamPort->Delete();

  return VTK_THREAD_RETURN_VALUE;
}


VTK_THREAD_RETURN_TYPE process_b( void *vtkNotUsed(arg) )
{
  vtkMultiProcessController *controller;
  int myid, otherid;
  
  //putenv("DISPLAY=:0.0");
  
  controller = vtkMultiProcessController::RegisterAndGetGlobalController(NULL);
  myid = controller->GetLocalProcessId();
  if (myid == 0)
    {
    otherid = 1;
    }
  else
    {
    otherid = 0;
    }

  vtkDownStreamPort *downStreamPort = vtkDownStreamPort::New();
  downStreamPort->SetUpStreamProcessId(otherid);
  downStreamPort->SetTag(999);

  vtkTexture *atext = vtkTexture::New();
  atext->SetInput(downStreamPort->GetImageDataOutput());
  //atext->SetInput(downStreamPort->GetStructuredPointsOutput());
  atext->InterpolateOn();

  vtkPlaneSource *plane = vtkPlaneSource::New();
  vtkPolyDataMapper  *mapper = vtkPolyDataMapper::New();
  mapper->SetInput(plane->GetOutput());
  
  vtkActor *actor = vtkActor::New();
  actor->SetMapper(mapper);
  actor->SetTexture(atext);
  
  // assign our actor to the renderer
  vtkRenderer *ren = vtkRenderer::New();
  ren->AddActor(actor);
  
  vtkRenderWindow *renWindow = vtkRenderWindow::New();
  renWindow->AddRenderer(ren);
  renWindow->SetSize( 300, 300 );

  vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
  iren->SetRenderWindow(renWindow);  
  
  // draw the resulting scene
  renWindow->Render();
  
  //  Begin mouse interaction
  iren->Start();
  
  // Clean up
  ren->Delete();
  renWindow->Delete();
  iren->Delete();
  downStreamPort->Delete();
  atext->Delete();
  plane->Delete();
  mapper->Delete();
  actor->Delete();

  return VTK_THREAD_RETURN_VALUE;
}


void main( int argc, char *argv[] )
{
  vtkMultiProcessController *controller;
  int myid;
  
  controller = vtkMultiProcessController::RegisterAndGetGlobalController(NULL);

  controller->Initialize(argc, argv);
  controller->SetNumberOfProcesses(2);
  controller->SetMultipleMethod(1, process_a, NULL);
  controller->SetMultipleMethod(0, process_b, NULL);
  controller->MultipleMethodExecute();

  controller->UnRegister(NULL);
}


