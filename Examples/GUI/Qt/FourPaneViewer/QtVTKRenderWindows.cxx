#include "ui_QtVTKRenderWindows.h"
#include "QtVTKRenderWindows.h"

#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include "vtkResliceImageViewer.h"
#include "vtkResliceCursorLineRepresentation.h"
#include "vtkResliceCursorThickLineRepresentation.h"
#include "vtkResliceCursorWidget.h"
#include "vtkResliceCursorActor.h"
#include "vtkResliceCursorPolyDataAlgorithm.h"
#include "vtkResliceCursor.h"
#include "vtkDICOMImageReader.h"
#include "vtkCellPicker.h"
#include "vtkProperty.h"
#include "vtkPlane.h"
#include "vtkImageData.h"
#include "vtkCommand.h"
#include "vtkPlaneSource.h"
#include "vtkLookupTable.h"
#include "vtkImageMapToWindowLevelColors.h"
#include "vtkInteractorStyleImage.h"
#include "vtkImageSlabReslice.h"


//----------------------------------------------------------------------------
class vtkResliceCursorCallback : public vtkCommand
{
public:
  static vtkResliceCursorCallback *New()
  { return new vtkResliceCursorCallback; }

  void Execute( vtkObject *caller, unsigned long ev,
                void *callData )
    {

    if (ev == vtkResliceCursorWidget::WindowLevelEvent ||
        ev == vtkCommand::WindowLevelEvent ||
        ev == vtkResliceCursorWidget::ResliceThicknessChangedEvent)
      {
      // Render everything
      for (int i = 0; i < 3; i++)
        {
        this->RCW[i]->Render();
        }
      this->IPW[0]->GetInteractor()->GetRenderWindow()->Render();
      return;
      }

    vtkImagePlaneWidget* ipw =
      dynamic_cast< vtkImagePlaneWidget* >( caller );
    if (ipw)
      {
      double* wl = static_cast<double*>( callData );

      if ( ipw == this->IPW[0] )
        {
        this->IPW[1]->SetWindowLevel(wl[0],wl[1],1);
        this->IPW[2]->SetWindowLevel(wl[0],wl[1],1);
        }
      else if( ipw == this->IPW[1] )
        {
        this->IPW[0]->SetWindowLevel(wl[0],wl[1],1);
        this->IPW[2]->SetWindowLevel(wl[0],wl[1],1);
        }
      else if (ipw == this->IPW[2])
        {
        this->IPW[0]->SetWindowLevel(wl[0],wl[1],1);
        this->IPW[1]->SetWindowLevel(wl[0],wl[1],1);
        }
      }

    vtkResliceCursorWidget *rcw = dynamic_cast<
      vtkResliceCursorWidget * >(caller);
    if (rcw)
      {
      vtkResliceCursorLineRepresentation *rep = dynamic_cast<
        vtkResliceCursorLineRepresentation * >(rcw->GetRepresentation());
      vtkResliceCursor *rc = rep->GetResliceCursorActor()->
                  GetCursorAlgorithm()->GetResliceCursor();
      for (int i = 0; i < 3; i++)
        {
        vtkPlaneSource *ps = static_cast< vtkPlaneSource * >(
            this->IPW[i]->GetPolyDataAlgorithm());
        ps->SetNormal(rc->GetPlane(i)->GetNormal());
        ps->SetCenter(rc->GetPlane(i)->GetOrigin());

        // If the reslice plane has modified, update it on the 3D widget
        this->IPW[i]->UpdatePlacement();
        }
      }

    // Render everything
    for (int i = 0; i < 3; i++)
      {
      this->RCW[i]->Render();
      }
    this->IPW[0]->GetInteractor()->GetRenderWindow()->Render();
    }

  vtkResliceCursorCallback() {}
  vtkImagePlaneWidget* IPW[3];
  vtkResliceCursorWidget *RCW[3];
};


QtVTKRenderWindows::QtVTKRenderWindows( int argc, char *argv[])
{
  this->ui = new Ui_QtVTKRenderWindows;
  this->ui->setupUi(this);

  vtkSmartPointer< vtkDICOMImageReader > reader =
    vtkSmartPointer< vtkDICOMImageReader >::New();
  reader->SetDirectoryName(argv[1]); //"/mnt/IOMegaFAT32USB/Kitware/VolView/src/VolViewData/Data/Osirix/AGECANONIX/Specials-1CoronaryCTA_with_spiral_CTA_pre/CorCTA-w-c-1.0-B20f");
  reader->Update();
  int imageDims[3];
  reader->GetOutput()->GetDimensions(imageDims);


  for (int i = 0; i < 3; i++)
    {
    riw[i] = vtkSmartPointer< vtkResliceImageViewer >::New();
    }

  this->ui->view1->SetRenderWindow(riw[0]->GetRenderWindow());
  riw[0]->SetupInteractor(
      this->ui->view1->GetRenderWindow()->GetInteractor());

  this->ui->view2->SetRenderWindow(riw[1]->GetRenderWindow());
  riw[1]->SetupInteractor(
      this->ui->view2->GetRenderWindow()->GetInteractor());

  this->ui->view3->SetRenderWindow(riw[2]->GetRenderWindow());
  riw[2]->SetupInteractor(
      this->ui->view3->GetRenderWindow()->GetInteractor());

  for (int i = 0; i < 3; i++)
    {      
    // make them all share the same reslice cursor object.
    vtkResliceCursorLineRepresentation *rep =
      vtkResliceCursorLineRepresentation::SafeDownCast(
          riw[i]->GetResliceCursorWidget()->GetRepresentation());
    rep->GetResliceCursorActor()->
      GetCursorAlgorithm()->SetResliceCursor(riw[0]->GetResliceCursor());

    rep->GetResliceCursorActor()->
      GetCursorAlgorithm()->SetReslicePlaneNormal(i);

    riw[i]->SetInput(reader->GetOutput());
    riw[i]->SetSliceOrientation(i);
    riw[i]->SetResliceModeToAxisAligned();
    }

  vtkSmartPointer<vtkCellPicker> picker =
    vtkSmartPointer<vtkCellPicker>::New();
  picker->SetTolerance(0.005);

  vtkSmartPointer<vtkProperty> ipwProp =
    vtkSmartPointer<vtkProperty>::New();

  vtkSmartPointer< vtkRenderer > ren =
    vtkSmartPointer< vtkRenderer >::New();

  this->ui->view4->GetRenderWindow()->AddRenderer(ren);
  vtkRenderWindowInteractor *iren = this->ui->view4->GetInteractor();

  for (int i = 0; i < 3; i++)
    {
    planeWidget[i] = vtkSmartPointer<vtkImagePlaneWidget>::New();
    planeWidget[i]->SetInteractor( iren );
    planeWidget[i]->SetPicker(picker);
    planeWidget[i]->RestrictPlaneToVolumeOn();
    double color[3] = {0, 0, 0};
    color[i] = 1;
    planeWidget[i]->GetPlaneProperty()->SetColor(color);

    color[0] /= 4.0;
    color[1] /= 4.0;
    color[2] /= 4.0;
    riw[i]->GetRenderer()->SetBackground( color );
    
    planeWidget[i]->SetTexturePlaneProperty(ipwProp);
    planeWidget[i]->TextureInterpolateOff();
    planeWidget[i]->SetResliceInterpolateToLinear();
    planeWidget[i]->SetInput(reader->GetOutput());
    planeWidget[i]->SetPlaneOrientation(i);
    planeWidget[i]->SetSliceIndex(imageDims[i]/2);
    planeWidget[i]->DisplayTextOn();
    planeWidget[i]->SetDefaultRenderer(ren);
    planeWidget[i]->SetWindowLevel(1358, -27);
    planeWidget[i]->On();
    planeWidget[i]->InteractionOn();
    }

  vtkSmartPointer<vtkResliceCursorCallback> cbk =
    vtkSmartPointer<vtkResliceCursorCallback>::New();

  for (int i = 0; i < 3; i++)
    {
    cbk->IPW[i] = planeWidget[i];
    cbk->RCW[i] = riw[i]->GetResliceCursorWidget();
    riw[i]->GetResliceCursorWidget()->AddObserver(
        vtkResliceCursorWidget::ResliceAxesChangedEvent, cbk );
    riw[i]->GetResliceCursorWidget()->AddObserver(
        vtkResliceCursorWidget::WindowLevelEvent, cbk );
    riw[i]->GetResliceCursorWidget()->AddObserver(
        vtkResliceCursorWidget::ResliceThicknessChangedEvent, cbk );
    riw[i]->GetResliceCursorWidget()->AddObserver(
        vtkResliceCursorWidget::ResetCursorEvent, cbk );
    riw[i]->GetInteractorStyle()->AddObserver(
        vtkCommand::WindowLevelEvent, cbk );

    // Make them all share the same color map.
    riw[i]->SetLookupTable(riw[0]->GetLookupTable());
    planeWidget[i]->GetColorMap()->SetLookupTable(riw[0]->GetLookupTable());
    }

  this->ui->view1->show();
  this->ui->view2->show();
  this->ui->view3->show();

  // Set up action signals and slots
  connect(this->ui->actionExit, SIGNAL(triggered()), this, SLOT(slotExit()));
  connect(this->ui->resliceModeCheckBox, SIGNAL(stateChanged(int)), this, SLOT(resliceMode(int)));
  connect(this->ui->thickModeCheckBox, SIGNAL(stateChanged(int)), this, SLOT(thickMode(int)));
  this->ui->thickModeCheckBox->setEnabled(0);

  connect(this->ui->radioButton_Max, SIGNAL(pressed()), this, SLOT(SetBlendModeToMaxIP()));
  connect(this->ui->radioButton_Min, SIGNAL(pressed()), this, SLOT(SetBlendModeToMinIP()));
  connect(this->ui->radioButton_Mean, SIGNAL(pressed()), this, SLOT(SetBlendModeToMeanIP()));
  this->ui->blendModeGroupBox->setEnabled(0);

  connect(this->ui->resetButton, SIGNAL(pressed()), this, SLOT(ResetViews()));
};

void QtVTKRenderWindows::slotExit()
{
  qApp->exit();
}

void QtVTKRenderWindows::resliceMode(int mode)
{
  this->ui->thickModeCheckBox->setEnabled(mode ? 1 : 0);
  this->ui->blendModeGroupBox->setEnabled(mode ? 1 : 0);

  for (int i = 0; i < 3; i++)
    {
    riw[i]->SetResliceMode(mode ? 1 : 0);
    riw[i]->GetRenderer()->ResetCamera();
    riw[i]->Render();
    }
}

void QtVTKRenderWindows::thickMode(int mode)
{
  for (int i = 0; i < 3; i++)
    {
    riw[i]->SetThickMode(mode ? 1 : 0);
    riw[i]->Render();
    }
}

void QtVTKRenderWindows::SetBlendMode(int m)
{
  for (int i = 0; i < 3; i++)
    {
    vtkImageSlabReslice *thickSlabReslice = vtkImageSlabReslice::SafeDownCast(
        vtkResliceCursorThickLineRepresentation::SafeDownCast(
          riw[i]->GetResliceCursorWidget()->GetRepresentation())->GetReslice());
    thickSlabReslice->SetBlendMode(m);
    riw[i]->Render();
    }
}

void QtVTKRenderWindows::SetBlendModeToMaxIP()
{
  this->SetBlendMode(VTK_IMAGESLAB_BLEND_MAX);
}

void QtVTKRenderWindows::SetBlendModeToMinIP()
{
  this->SetBlendMode(VTK_IMAGESLAB_BLEND_MIN);
}

void QtVTKRenderWindows::SetBlendModeToMeanIP()
{
  this->SetBlendMode(VTK_IMAGESLAB_BLEND_MEAN);
}

void QtVTKRenderWindows::ResetViews()
{
  // Reset the reslice image views
  for (int i = 0; i < 3; i++)
    {
    riw[i]->Reset();
    }

  // Also sync the Image plane widget on the 3D top right view with any
  // changes to the reslice cursor.
  for (int i = 0; i < 3; i++)
    {
    vtkPlaneSource *ps = static_cast< vtkPlaneSource * >(
        planeWidget[i]->GetPolyDataAlgorithm());
    ps->SetNormal(riw[0]->GetResliceCursor()->GetPlane(i)->GetNormal());
    ps->SetCenter(riw[0]->GetResliceCursor()->GetPlane(i)->GetOrigin());

    // If the reslice plane has modified, update it on the 3D widget
    this->planeWidget[i]->UpdatePlacement();
    }

  // Render in response to changes.
  this->Render();
}

void QtVTKRenderWindows::Render()
{
  for (int i = 0; i < 3; i++)
    {
    riw[i]->Render();
    }
  this->ui->view3->GetRenderWindow()->Render();
}
