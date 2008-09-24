/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkQtChartColorStyleGenerator.h

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

/// \file vtkQtChartColorStyleGenerator.h
/// \date September 22, 2008

#ifndef _vtkQtChartColorStyleGenerator_h
#define _vtkQtChartColorStyleGenerator_h

#include "vtkQtChartExport.h"
#include "vtkQtChartStyleGenerator.h"
#include "vtkQtChartColors.h" // Needed for enum.

class vtkQtChartColorStyleGeneratorInternal;


/// \class vtkQtChartColorStyleGenerator
/// \brief
///   The vtkQtChartColorStyleGenerator class generates series drawing
///   options using color and pen style.
class VTKQTCHART_EXPORT vtkQtChartColorStyleGenerator :
  public vtkQtChartStyleGenerator
{
  Q_OBJECT

public:
  /// \brief
  ///   Creates a color/style generator.
  /// \param parent The parent object.
  /// \param scheme The initial color scheme.
  vtkQtChartColorStyleGenerator(QObject *parent=0,
    vtkQtChartColors::ColorScheme scheme=vtkQtChartColors::Spectrum);
  virtual ~vtkQtChartColorStyleGenerator();

  /// \name vtkQtChartStyleGenerator Methods
  //@{
  /// \brief
  ///   Gets the brush for the specified series index.
  ///
  /// If the index is greater than the internal color list, the index
  /// will be wrapped around repeating the colors.
  ///
  /// \param index The series index.
  /// \return
  ///   The brush for the specified series index.
  virtual QBrush getSeriesBrush(int index) const;

  /// \brief
  ///   Gets the pen for the specified series index.
  ///
  /// If the index is greater than the internal color list, the index
  /// will be wrapped around repeating the colors. The repeated
  /// colors will have the next pen style in the list.
  ///
  /// \param index The series index.
  /// \return
  ///   The pen for the specified series index.
  virtual QPen getSeriesPen(int index) const;
  //@}

  /// \name Color Methods
  //@{
  /// \brief
  ///   Gets the list of colors.
  /// \return
  ///   A pointer to the list of colors.
  vtkQtChartColors *getColors() {return this->Colors;}

  /// \brief
  ///   Gets the list of colors.
  /// \return
  ///   A pointer to the list of colors.
  const vtkQtChartColors *getColors() const {return this->Colors;}
  //@}

  /// \name Pen Style Methods
  //@{
  /// \brief
  ///   Gets the number of pen styles in the style list.
  /// \return
  ///   The number of pen styles in the style list.
  int getNumberOfStyles() const;

  /// \brief
  ///   Gets a pen style from the pen styles list.
  ///
  /// This method provides access to the list of styles. If the index
  /// is out of range, a default style will be returned.
  ///
  /// \param index The list index for the style.
  /// \return
  ///   The pen style for the given index.
  /// \sa vtkQtChartStyleGenerator::getSeriesPen(int)
  Qt::PenStyle getPenStyle(int index) const;

  /// \brief
  ///   Sets the pen style for the given index.
  ///
  /// This method does nothing if the index is out of range.
  ///
  /// \param index Which pen style to modify.
  /// \param style The new pen style.
  void setPenStyle(int index, Qt::PenStyle style);

  /// Clears the list of pen styles.
  void clearPenStyles();

  /// \brief
  ///   Adds a pen style to the list of pen styles.
  /// \param style The new pen style to add.
  void addPenStyle(Qt::PenStyle style);

  /// \brief
  ///   Inserts a new pen style into the list of pen styles.
  /// \param index Where to insert the new pen style.
  /// \param style The new pen style to insert.
  void insertPenStyle(int index, Qt::PenStyle style);

  /// \brief
  ///   Removes the pen style for the given index.
  /// \param index Which pen style to remove from the list.
  void removePenStyle(int index);
  //@}

private:
  /// Stores the list of pen styles.
  vtkQtChartColorStyleGeneratorInternal *Internal;
  vtkQtChartColors *Colors; ///< Stores the color list.
};

#endif
