/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkQtTreeView.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
/*-------------------------------------------------------------------------
  Copyright 2008 Sandia Corporation.
  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
  the U.S. Government retains certain rights in this software.
-------------------------------------------------------------------------*/
// .NAME vtkQtTreeView - A VTK view based on a Qt tree view.
//
// .SECTION Description
// vtkQtTreeView is a VTK view using an underlying QTreeView. 
//
// .SECTION Thanks
// Thanks to Brian Wylie from Sandia National Laboratories for implementing
// this class

#ifndef __vtkQtTreeView_h
#define __vtkQtTreeView_h

#include "QVTKWin32Header.h"
#include "vtkQtView.h"

#include <QPointer>
#include "vtkQtAbstractModelAdapter.h"

class QItemSelection;
class QModelIndex;
class QTreeView;
class vtkQtTreeModelAdapter;

class QVTK_EXPORT vtkQtTreeView : public vtkQtView
{
Q_OBJECT

signals:
  void expanded(const QModelIndex&);
  void collapsed(const QModelIndex&);

public:
  static vtkQtTreeView *New();
  vtkTypeRevisionMacro(vtkQtTreeView, vtkQtView);
  void PrintSelf(ostream& os, vtkIndent indent);
  
  // Description:
  // Get the main container of this view (a  QWidget).
  // The application typically places the view with a call
  // to GetWidget(): something like this
  // this->ui->box->layout()->addWidget(this->View->GetWidget());
  virtual QWidget* GetWidget();

  // Description:
  // Pointer to the internal model adapter used convert the
  // vtkDataObject to a QAbstractItemModel.
  vtkQtAbstractModelAdapter* GetItemModelAdapter();

  void ExpandAll();

  // Description:
  // Updates the view.
  virtual void Update();

protected:
  vtkQtTreeView();
  ~vtkQtTreeView();

  // Description:
  // Connects the algorithm output to the internal pipeline.
  // This view only supports a single representation.
  virtual void AddInputConnection( int port, int index,
    vtkAlgorithmOutput* conn,
    vtkAlgorithmOutput* selectionConn);
  
  // Description:
  // Removes the algorithm output from the internal pipeline.
  virtual void RemoveInputConnection( int port, int index,
    vtkAlgorithmOutput* conn,
    vtkAlgorithmOutput* selectionConn);

  // Description:
  // We need to keep track of whether were in selection mode
  bool Selecting;
  
  QPointer<QTreeView> TreeView;
  vtkQtTreeModelAdapter* TreeAdapter;

private slots:
  void slotSelectionChanged(const QItemSelection&,const QItemSelection&);

private:
  vtkQtTreeView(const vtkQtTreeView&);  // Not implemented.
  void operator=(const vtkQtTreeView&);  // Not implemented.
  
};

#endif
