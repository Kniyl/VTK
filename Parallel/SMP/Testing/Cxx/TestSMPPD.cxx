#include "vtkActor.h"
#include "vtkAssignAttribute.h"
#include "vtkContourFilter.h"
#include "vtkDataSetMapper.h"
#include "vtkElevationFilter.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRTAnalyticSource.h"
#include "vtkSMPContourFilter.h"
#include "vtkSMPMergePoints.h"
#include "vtkSMPMinMaxTree.h"
#include "vtkSMPTransform.h"
#include "vtkTestUtilities.h"
#include "vtkThreshold.h"
#include "vtkTimerLog.h"
#include "vtkTransform.h"
#include "vtkTransformFilter.h"
#include "vtkUnstructuredGrid.h"

#include <cstdlib>

const int NumRepetitions = 1;
const int WaveSize = 100;

void testPD(vtkContourFilter *isosurface)
{
  vtkTimerLog *timer = vtkTimerLog::New();
  isosurface->SetInputArrayToProcess(0,0,0, vtkDataObject::FIELD_ASSOCIATION_POINTS, "Elevation");
  isosurface->GenerateValues( 150, 0.0, 1.0 );
  isosurface->UseScalarTreeOn();
  cerr << "update " << isosurface->GetClassName() << endl;
  double t, t0,t1;
  t = 0;
  for (int i = 0; i < NumRepetitions; ++i)
    {
    isosurface->Modified();
    t0 = timer->GetUniversalTime();
    isosurface->Update();
    t1 = timer->GetUniversalTime();
    t += t1-t0;
    }
  cerr << (t)/NumRepetitions << endl;
  timer->Delete();
}

int TestSMPPD( int argc, char * argv [] )
{
  vtkTimerLog *timer = vtkTimerLog::New();
  double t0, t1;

  int threads = 0;
  if (argc > 1)
    {
    threads = atoi(argv[1]);
    }
  bool sequential = (threads==0);

  vtkRTAnalyticSource* wavelet = vtkRTAnalyticSource::New();
  wavelet->SetWholeExtent(-WaveSize,WaveSize,-WaveSize,WaveSize,-WaveSize,WaveSize);
  wavelet->SetCenter(0,0,0);
  wavelet->SetMaximum(255);
  wavelet->SetXFreq(60);
  wavelet->SetYFreq(30);
  wavelet->SetZFreq(40);
  wavelet->SetXMag(10);
  wavelet->SetYMag(18);
  wavelet->SetZMag(5);
  wavelet->SetStandardDeviation(0.5);
  wavelet->SetSubsampleRate(1);

  vtkContourFilter *cf = vtkContourFilter::New();
  cf->SetInputConnection(wavelet->GetOutputPort());
  cf->SetNumberOfContours(1);
  cf->SetValue(0, 150);
  wavelet->Delete();

  double b[6];
  cf->Update();
  cf->GetOutput()->GetBounds(b);

  vtkElevationFilter *ef = vtkElevationFilter::New();
  ef->SetInputConnection(cf->GetOutputPort());
  ef->SetLowPoint(b[0], b[2], b[4]);
  ef->SetHighPoint(b[1], b[3], b[5]);

  vtkAssignAttribute* aa = vtkAssignAttribute::New();
  aa->SetInputConnection(ef->GetOutputPort());
  aa->Assign("Elevation", vtkDataSetAttributes::SCALARS, vtkAssignAttribute::POINT_DATA);
  ef->Delete();

  vtkTransformFilter* transform = vtkTransformFilter::New();
  transform->SetInputConnection( aa->GetOutputPort() );
  aa->Delete();

  if ( sequential )
    {
    vtkTransform* t = vtkTransform::New();
    t->Scale( 1, 2, 3 );
    transform->SetTransform( t );
    cerr << "update " << t->GetClassName() << endl;
    t->Delete();
    }
  else
    {
    vtkSMPTransform* t = vtkSMPTransform::New();
    t->Scale( 1, 2, 3 );
    transform->SetTransform( t );
    cerr << "update " << t->GetClassName() << endl;
    t->Delete();
    }

  t0 = timer->GetUniversalTime();
  for (int i = 0; i < NumRepetitions; ++i)
    {
    transform->Modified();
    transform->Update();
    }
  t1 = timer->GetUniversalTime();
  cerr << (t1-t0)/NumRepetitions << endl;

  /* === Testing contour filter === */

  vtkContourFilter* isosurface1 = vtkContourFilter::New();
  isosurface1->SetInputConnection( transform->GetOutputPort() );

  vtkContourFilter* isosurface2 = vtkSMPContourFilter::New();
  vtkSMPMergePoints* locator = vtkSMPMergePoints::New();
  isosurface2->SetLocator( locator );
  locator->Delete();
  vtkSMPMinMaxTree* tree = vtkSMPMinMaxTree::New();
  isosurface2->SetScalarTree(tree);
  tree->Delete();
  isosurface2->SetInputConnection( transform->GetOutputPort() );

  testPD(isosurface1);
  testPD(isosurface2);

  /* === Watching outputs === */

  isosurface1->GetOutput()->GetBounds(b);

  vtkDataSetMapper* map1 = vtkDataSetMapper::New();
  map1->SetInputConnection( isosurface1->GetOutputPort() );
  vtkActor* actor1 = vtkActor::New();
  actor1->SetMapper( map1 );
  map1->Delete();

  vtkDataSetMapper* map2 = vtkDataSetMapper::New();
  map2->SetInputConnection( isosurface2->GetOutputPort() );
  vtkActor* actor2 = vtkActor::New();
  actor2->SetMapper( map2 );
  map2->Delete();
  actor2->AddPosition( (b[1]-b[0])*1.1, 0., 0. );

  vtkRenderer* viewport = vtkRenderer::New();
  viewport->SetBackground( .5, .5, .5 );
  viewport->AddActor( actor1 );
  viewport->AddActor( actor2 );
  actor1->Delete();
  actor2->Delete();

  vtkRenderWindow* window = vtkRenderWindow::New();
  window->AddRenderer( viewport );
  viewport->Delete();

  vtkRenderWindowInteractor* iren = vtkRenderWindowInteractor::New();
  iren->SetRenderWindow( window );
  window->Delete();

  iren->Initialize();

  iren->Start();

  iren->Delete();

  isosurface1->Delete();
  isosurface2->Delete();

  timer->Delete();
  return 0;
}
