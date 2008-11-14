/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkQtStackedChart.cxx

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

/// \file vtkQtStackedChart.cxx
/// \date February 27, 2008

#ifdef _MSC_VER
// Disable warnings that Qt headers give.
#pragma warning(disable:4127)
#endif

#include "vtkQtStackedChart.h"

#include "vtkQtChartAxis.h"
#include "vtkQtChartAxisCornerDomain.h"
#include "vtkQtChartAxisDomain.h"
#include "vtkQtChartAxisLayer.h"
#include "vtkQtChartAxisOptions.h"
#include "vtkQtChartColors.h"
#include "vtkQtChartContentsArea.h"
#include "vtkQtChartContentsSpace.h"
#include "vtkQtChartHelpFormatter.h"
#include "vtkQtChartLayerDomain.h"
#include "vtkQtChartQuad.h"
#include "vtkQtChartShapeLocator.h"
#include "vtkQtChartSeriesDomain.h"
#include "vtkQtChartSeriesDomainGroup.h"
#include "vtkQtChartSeriesModel.h"
#include "vtkQtChartSeriesSelection.h"
#include "vtkQtChartSeriesSelectionModel.h"
#include "vtkQtChartArea.h"
#include "vtkQtStackedChartOptions.h"
#include "vtkQtStackedChartSeriesOptions.h"

#include <QBrush>
#include <QGraphicsScene>
#include <QList>
#include <QPen>
#include <QPointF>
#include <QPolygonF>
#include <QStyleOptionGraphicsItem>
#include <QVector>


class vtkQtStackedChartSeries
{
public:
  vtkQtStackedChartSeries(QPolygonF *polygon=0);
  vtkQtStackedChartSeries(const vtkQtStackedChartSeries &other);
  ~vtkQtStackedChartSeries();

  vtkQtStackedChartSeries &operator=(const vtkQtStackedChartSeries &other);

  void setMapping(int group, int index);

  void updateGradient();

  void clearQuads();
  void clearHighlights();

public:
  QPolygonF *Polygon;
  QList<vtkQtChartQuad *> Quads;
  QList<QPolygonF *> Highlights;
  QPointF Gradient1;
  QPointF Gradient2;
  int Group;
  int Index;
  bool IsHighlighted;
};


class vtkQtStackedChartSeriesGroup
{
public:
  vtkQtStackedChartSeriesGroup();
  vtkQtStackedChartSeriesGroup(const vtkQtStackedChartSeriesGroup &other);
  ~vtkQtStackedChartSeriesGroup() {}

  vtkQtStackedChartSeriesGroup &operator=(
      const vtkQtStackedChartSeriesGroup &other);

public:
  QVector<QVector<double> > Data;
  QList<QList<vtkQtChartShape *> > Shapes;
};


class vtkQtStackedChartDomainGroup : public vtkQtChartSeriesDomainGroup
{
public:
  vtkQtStackedChartDomainGroup();
  virtual ~vtkQtStackedChartDomainGroup();

  virtual void clear();

protected:
  virtual void insertGroup(int group);
  virtual void removeGroup(int group);

private:
  void cleanUp();

public:
  QList<vtkQtStackedChartSeriesGroup *> Tables;
};


class vtkQtStackedChartInternal
{
public:
  vtkQtStackedChartInternal();
  ~vtkQtStackedChartInternal();

  QPointF getMidPoint(const QPointF &point1, const QPointF &point2) const;

  QList<vtkQtStackedChartSeries *> Series;
  vtkQtChartAxisCornerDomain Domain;
  vtkQtStackedChartDomainGroup Groups;
  vtkQtChartShapeLocator QuadTree;
  QRectF Bounds;
  int CurrentGroup;
};


//-----------------------------------------------------------------------------
vtkQtStackedChartSeries::vtkQtStackedChartSeries(QPolygonF *polygon)
  : Quads(), Highlights(), Gradient1(), Gradient2()
{
  this->Polygon = polygon;
  this->Group = -1;
  this->Index = -1;
  this->IsHighlighted = false;
}

vtkQtStackedChartSeries::vtkQtStackedChartSeries(
    const vtkQtStackedChartSeries &other)
  : Highlights(), Gradient1(other.Gradient1), Gradient2(other.Gradient2)
{
  this->Polygon = 0;
  this->Group = other.Group;
  this->Index = other.Index;
  this->IsHighlighted = other.IsHighlighted;
  if(other.Polygon)
    {
    this->Polygon = new QPolygonF(*other.Polygon);
    }

  // TODO: Copy the quads and highlights.
}

vtkQtStackedChartSeries::~vtkQtStackedChartSeries()
{
  this->clearQuads();
  this->clearHighlights();
  if(this->Polygon)
    {
    delete this->Polygon;
    }
}

vtkQtStackedChartSeries &vtkQtStackedChartSeries::operator=(
    const vtkQtStackedChartSeries &other)
{
  this->Gradient1 = other.Gradient1;
  this->Gradient2 = other.Gradient2;
  this->Group = other.Group;
  this->Index = other.Index;
  this->IsHighlighted = other.IsHighlighted;
  if(this->Polygon)
    {
    if(other.Polygon)
      {
      *this->Polygon = *other.Polygon;
      }
    else
      {
      delete this->Polygon;
      this->Polygon = 0;
      }
    }
  else if(other.Polygon)
    {
    this->Polygon = new QPolygonF(*other.Polygon);
    }

  // TODO: Copy the quads and highlights.
  this->clearQuads();
  this->clearHighlights();

  return *this;
}

void vtkQtStackedChartSeries::setMapping(int group, int index)
{
  this->Group = group;
  this->Index = index;
}

void vtkQtStackedChartSeries::updateGradient()
{
  QRectF bounds = this->Polygon->boundingRect();
  float center = bounds.center().x();
  this->Gradient1.setX(center);
  this->Gradient1.setY(bounds.top());
  this->Gradient2.setX(center);
  this->Gradient2.setY(bounds.bottom());
}

void vtkQtStackedChartSeries::clearHighlights()
{
  QList<QPolygonF *>::Iterator iter = this->Highlights.begin();
  for( ; iter != this->Highlights.end(); ++iter)
    {
    delete *iter;
    }

  this->Highlights.clear();
}

void vtkQtStackedChartSeries::clearQuads()
{
  QList<vtkQtChartQuad *>::Iterator iter = this->Quads.begin();
  for( ; iter != this->Quads.end(); ++iter)
    {
    delete *iter;
    }

  this->Quads.clear();
}


//-----------------------------------------------------------------------------
vtkQtStackedChartSeriesGroup::vtkQtStackedChartSeriesGroup()
  : Data(), Shapes()
{
}

vtkQtStackedChartSeriesGroup::vtkQtStackedChartSeriesGroup(
    const vtkQtStackedChartSeriesGroup &other)
  : Data(other.Data), Shapes(other.Shapes)
{
}

vtkQtStackedChartSeriesGroup &vtkQtStackedChartSeriesGroup::operator=(
    const vtkQtStackedChartSeriesGroup &other)
{
  this->Data = other.Data;
  this->Shapes = other.Shapes;
  return *this;
}


//-----------------------------------------------------------------------------
vtkQtStackedChartDomainGroup::vtkQtStackedChartDomainGroup()
  : vtkQtChartSeriesDomainGroup(true), Tables()
{
}

vtkQtStackedChartDomainGroup::~vtkQtStackedChartDomainGroup()
{
  this->cleanUp();
}

void vtkQtStackedChartDomainGroup::clear()
{
  vtkQtChartSeriesDomainGroup::clear();
  this->cleanUp();
  this->Tables.clear();
}

void vtkQtStackedChartDomainGroup::insertGroup(int group)
{
  vtkQtChartSeriesDomainGroup::insertGroup(group);
  this->Tables.insert(group, new vtkQtStackedChartSeriesGroup());
}

void vtkQtStackedChartDomainGroup::removeGroup(int group)
{
  vtkQtChartSeriesDomainGroup::removeGroup(group);
  delete this->Tables.takeAt(group);
}

void vtkQtStackedChartDomainGroup::cleanUp()
{
  QList<vtkQtStackedChartSeriesGroup *>::Iterator iter = this->Tables.begin();
  for( ; iter != this->Tables.end(); ++iter)
    {
    delete *iter;
    }
}


//-----------------------------------------------------------------------------
vtkQtStackedChartInternal::vtkQtStackedChartInternal()
  : Series(), Domain(), Groups(), QuadTree(), Bounds()
{
  this->CurrentGroup = -1;

  this->Domain.setVerticalPreferences(false, true, false);
}

vtkQtStackedChartInternal::~vtkQtStackedChartInternal()
{
  // Clean up the remaining series items.
  QList<vtkQtStackedChartSeries *>::Iterator iter = this->Series.begin();
  for( ; iter != this->Series.end(); ++iter)
    {
    delete *iter;
    }
}

QPointF vtkQtStackedChartInternal::getMidPoint(const QPointF &point1,
    const QPointF &point2) const
{
  return QPointF(
      (point2.x() + point1.x()) * 0.5, (point2.y() + point1.y()) * 0.5);
}


//-----------------------------------------------------------------------------
vtkQtStackedChart::vtkQtStackedChart()
  : vtkQtChartSeriesLayer(false)
{
  this->Internal = new vtkQtStackedChartInternal();
  this->Options = new vtkQtStackedChartOptions(this);
  this->InModelChange = false;
  this->BuildNeeded = false;

  // Listen for option changes.
  this->connect(this->Options, SIGNAL(axesCornerChanged()),
      this, SLOT(handleAxesCornerChange()));
  this->connect(this->Options, SIGNAL(sumationChanged()),
      this, SLOT(handleSumationChange()));
  this->connect(this->Options, SIGNAL(gradientChanged()),
      this, SLOT(handleGradientChange()));

  // Listen for selection changes.
  this->connect(this->Selection,
      SIGNAL(selectionChanged(const vtkQtChartSeriesSelection &)),
      this, SLOT(updateHighlights()));
}

vtkQtStackedChart::~vtkQtStackedChart()
{
  delete this->Internal;
}

void vtkQtStackedChart::setChartArea(vtkQtChartArea *area)
{
  vtkQtChartSeriesLayer::setChartArea(area);
  this->reset();
}

void vtkQtStackedChart::setModel(vtkQtChartSeriesModel *model)
{
  if(this->Model)
    {
    // Disconnect from the previous model's signals.
    this->disconnect(this->Model, 0, this, 0);
    }

  vtkQtChartSeriesLayer::setModel(model);
  if(this->Model)
    {
    // Listen for model changes.
    this->connect(this->Model, SIGNAL(modelReset()), this, SLOT(reset()));
    this->connect(this->Model, SIGNAL(seriesAboutToBeInserted(int, int)),
        this, SLOT(prepareSeriesInsert(int, int)));
    this->connect(this->Model, SIGNAL(seriesInserted(int, int)),
        this, SLOT(insertSeries(int, int)));
    this->connect(this->Model, SIGNAL(seriesAboutToBeRemoved(int, int)),
        this, SLOT(startSeriesRemoval(int, int)));
    this->connect(this->Model, SIGNAL(seriesRemoved(int, int)),
        this, SLOT(finishSeriesRemoval(int, int)));
    }

  // Reset the view items for the new model.
  this->reset();
}

void vtkQtStackedChart::setOptions(const vtkQtStackedChartOptions &options)
{
  // Copy the new options. The chart will collapse the layout signals.
  this->Options->setSumNormalized(options.isSumNormalized());
  this->Options->setGradientDisplayed(options.isGradientDislpayed());
  this->Options->setAxesCorner(options.getAxesCorner());
  this->Options->getHelpFormat()->setFormat(
      options.getHelpFormat()->getFormat());
}

vtkQtStackedChartSeriesOptions *vtkQtStackedChart::getStackedSeriesOptions(
    int series) const
{
  return qobject_cast<vtkQtStackedChartSeriesOptions *>(
      this->getSeriesOptions(series));
}

void vtkQtStackedChart::getLayerDomain(vtkQtChartLayerDomain &domain) const
{
  domain.mergeDomain(this->Internal->Domain, this->Options->getAxesCorner());
}

void vtkQtStackedChart::layoutChart(const QRectF &area)
{
  // Update the position and bounds.
  this->prepareGeometryChange();
  this->Internal->Bounds.setSize(area.size());
  this->setPos(area.topLeft());
  if(this->Internal->Series.size() == 0)
    {
    return;
    }

  // Get the axis layer to get the axes and domains.
  vtkQtChartAxisLayer *layer = this->ChartArea->getAxisLayer();
  vtkQtChartAxis *xAxis = layer->getHorizontalAxis(
      this->Options->getAxesCorner());
  vtkQtChartAxis *yAxis = layer->getVerticalAxis(
      this->Options->getAxesCorner());

  int seriesGroup;
  const vtkQtChartSeriesDomain *seriesDomain =
      this->Internal->Domain.getDomain(xAxis->getAxisDomain(),
      yAxis->getAxisDomain(), &seriesGroup);

  // Get the x-axis minimum and maximum locations.
  float zero = yAxis->getZeroPixel();
  bool isRange = false;
  QList<QVariant> xDomain;
  QList<int> seriesList;
  vtkQtStackedChartSeriesGroup *group = 0;
  if(seriesDomain)
    {
    seriesList = this->Internal->Groups.getGroup(seriesGroup);
    xDomain = seriesDomain->getXDomain().getDomain(isRange);
    group = this->Internal->Groups.Tables[seriesGroup];
    }

  int i = 0;
  QList<vtkQtStackedChartSeries *>::Iterator iter =
      this->Internal->Series.begin();
  for(int series = 0; iter != this->Internal->Series.end(); ++iter, ++series)
    {
    QPolygonF *polygon = (*iter)->Polygon;
    if(polygon && seriesList.contains(series))
      {
      polygon->clear();
      int j = 0;
      int half = group->Data[i].size();
      for( ; j < half; j++)
        {
        polygon->append(QPointF(xAxis->getPixel(xDomain[j]), yAxis->getPixel(
            QVariant(group->Data[i][j]))));
        }

      j = half - 1;
      if(i == 0)
        {
        for( ; j >= 0; j--)
          {
          polygon->append(QPointF(xAxis->getPixel(xDomain[j]), zero));
          }
        }
      else
        {
        int k = i - 1;
        for( ; j >= 0; j--)
          {
          polygon->append(QPointF(xAxis->getPixel(xDomain[j]), yAxis->getPixel(
              QVariant(group->Data[k][j]))));
          }
        }

      // Build the series quads from the polygon outline.
      int total = polygon->size();
      QList<vtkQtChartQuad *>::Iterator jter = (*iter)->Quads.begin();
      for(j = 1; j < half && jter != (*iter)->Quads.end(); ++j, ++jter)
        {
        vtkQtChartQuad *left = *jter;
        ++jter;
        if(jter == (*iter)->Quads.end())
          {
          break;
          }

        // Find the midpoints of the line segments.
        vtkQtChartQuad *right = *jter;
        QPointF midTop = this->Internal->getMidPoint(
            (*polygon)[j - 1], (*polygon)[j]);
        QPointF midBottom = this->Internal->getMidPoint(
            (*polygon)[total - j], (*polygon)[total - j - 1]);

        // Set the quad points.
        left->setPoint(0, (*polygon)[j - 1]);
        left->setPoint(1, midTop);
        left->setPoint(2, midBottom);
        left->setPoint(3, (*polygon)[total - j]);

        right->setPoint(0, midTop);
        right->setPoint(1, (*polygon)[j]);
        right->setPoint(2, (*polygon)[total - j - 1]);
        right->setPoint(3, midBottom);
        }

      // Increment the data table index for the next series.
      i++;

      // Set up the series gradient if needed.
      if(this->Options->isGradientDislpayed())
        {
        (*iter)->updateGradient();
        }
      }
    }

  // Layout the highlights.
  this->layoutHighlights();

  // Update the quad tree.
  if(seriesDomain)
    {
    if(this->ChartArea->isInteractivelyResizing())
      {
      this->BuildNeeded = true;
      }
    else
      {
      this->buildQuadTree(seriesGroup);
      }
    }
}

bool vtkQtStackedChart::getHelpText(const QPointF &point, QString &text)
{
  vtkQtChartSeriesSelection selection;
  this->getPointsAt(point, selection);
  if(!selection.isEmpty())
    {
    // Use the axis options to format the data.
    vtkQtChartAxisLayer *layer = this->ChartArea->getAxisLayer();
    vtkQtChartAxisOptions *xAxis = layer->getHorizontalAxis(
        this->Options->getAxesCorner())->getOptions();
    vtkQtChartAxisOptions *yAxis = layer->getVerticalAxis(
        this->Options->getAxesCorner())->getOptions();

    // Use the x-axis domain and the table for the data values.
    const QList<vtkQtChartSeriesSelectionItem> &points = selection.getPoints();
    int series = points[0].Series;
    vtkQtStackedChartSeries *item = this->Internal->Series[series];
    const vtkQtChartSeriesDomain *seriesDomain =
        this->Internal->Domain.getDomain(item->Group);
    bool isRange = false;
    int index = points[0].Points[0].first;
    QStringList args;
    args.append(xAxis->formatValue(
        seriesDomain->getXDomain().getDomain(isRange)[index]));
    vtkQtStackedChartSeriesGroup *group =
        this->Internal->Groups.Tables[item->Group];
    args.append(yAxis->formatValue(QVariant(
        group->Data[item->Index][index])));
    if(item->Index > 0)
      {
      double value = group->Data[item->Index][index] -
          group->Data[item->Index - 1][index];
      args.append(yAxis->formatValue(QVariant(value)));
      }
    else
      {
      args.append(args[1]);
      }

    text = this->Options->getHelpFormat()->getHelpText(
        this->Model->getSeriesName(series).toString(), args);
    return true;
    }

  return false;
}

void vtkQtStackedChart::finishInteractiveResize()
{
  if(this->BuildNeeded)
    {
    // Get the axis layer to get the axes and domains.
    vtkQtChartAxisLayer *layer = this->ChartArea->getAxisLayer();
    vtkQtChartAxis *xAxis = layer->getHorizontalAxis(
        this->Options->getAxesCorner());
    vtkQtChartAxis *yAxis = layer->getVerticalAxis(
        this->Options->getAxesCorner());

    int seriesGroup;
    const vtkQtChartSeriesDomain *seriesDomain =
        this->Internal->Domain.getDomain(xAxis->getAxisDomain(),
        yAxis->getAxisDomain(), &seriesGroup);
    if(seriesDomain)
      {
      this->buildQuadTree(seriesGroup);
      }
    }
}

void vtkQtStackedChart::getSeriesAt(const QPointF &point,
    vtkQtChartSeriesSelection &selection) const
{
  // Translate the point to contents coordinates.
  QPointF local = point;
  this->ChartArea->getContentsSpace()->translateToLayerContents(local);

  // Get the selected quads from the tree.
  vtkQtChartIndexRangeList indexes;
  QList<vtkQtChartShape *> shapes = this->Internal->QuadTree.getItemsAt(local);
  if(shapes.size() > 0)
    {
    // Add the series to the selection.
    int series = shapes.first()->getSeries();
    indexes.append(vtkQtChartIndexRange(series, series));
    }

  selection.setSeries(indexes);
}

void vtkQtStackedChart::getPointsAt(const QPointF &point,
    vtkQtChartSeriesSelection &selection) const
{
  // Translate the point to contents coordinates.
  QPointF local = point;
  this->ChartArea->getContentsSpace()->translateToLayerContents(local);

  // Get the bar index from the search tree.
  QList<vtkQtChartSeriesSelectionItem> indexes;
  QList<vtkQtChartShape *> shapes = this->Internal->QuadTree.getItemsAt(local);
  if(shapes.size() > 0)
    {
    vtkQtChartSeriesSelectionItem item(shapes.first()->getSeries());
    int index = shapes.first()->getIndex();
    item.Points.append(vtkQtChartIndexRange(index, index));
    indexes.append(item);
    }

  selection.setPoints(indexes);
}

void vtkQtStackedChart::getSeriesIn(const QRectF &area,
    vtkQtChartSeriesSelection &selection) const
{
  // Translate the rectangle to contents coordinates.
  QRectF local = area;
  this->ChartArea->getContentsSpace()->translateToLayerContents(local);

  // Get the list of bar indexes from the bar tree.
  vtkQtChartIndexRangeList indexes;
  QList<vtkQtChartShape *> shapes = this->Internal->QuadTree.getItemsIn(local);
  QList<vtkQtChartShape *>::Iterator shape = shapes.begin();
  for( ; shape != shapes.end(); ++shape)
    {
    // Add the series to the selection.
    int series = (*shape)->getSeries();
    indexes.append(vtkQtChartIndexRange(series, series));
    }

  selection.setSeries(indexes);
}

void vtkQtStackedChart::getPointsIn(const QRectF &area,
    vtkQtChartSeriesSelection &selection) const
{
  // Translate the rectangle to contents coordinates.
  QRectF local = area;
  this->ChartArea->getContentsSpace()->translateToLayerContents(local);

  // Get the list of bar indexes from the bar tree.
  QList<vtkQtChartSeriesSelectionItem> indexes;
  QList<vtkQtChartShape *> shapes = this->Internal->QuadTree.getItemsIn(local);
  QList<vtkQtChartShape *>::Iterator shape = shapes.begin();
  for( ; shape != shapes.end(); ++shape)
    {
    vtkQtChartSeriesSelectionItem item((*shape)->getSeries());
    int index = (*shape)->getIndex();
    item.Points.append(vtkQtChartIndexRange(index, index));
    indexes.append(item);
    }

  selection.setPoints(indexes);
}

QRectF vtkQtStackedChart::boundingRect() const
{
  return this->Internal->Bounds;
}

void vtkQtStackedChart::paint(QPainter *painter,
    const QStyleOptionGraphicsItem *option, QWidget *)
{
  if(!this->ChartArea)
    {
    return;
    }

  // Use the exposed rectangle from the option object to determine
  // which series to draw.
  vtkQtChartContentsSpace *space = this->ChartArea->getContentsSpace();
  QRectF area = option->exposedRect.translated(space->getXOffset(),
      space->getYOffset());

  // Get the axis layer to get the axes and domain priority.
  vtkQtChartAxisLayer *layer = this->ChartArea->getAxisLayer();
  vtkQtChartLayer::AxesCorner corner = this->Options->getAxesCorner();
  vtkQtChartAxis *xAxis = layer->getHorizontalAxis(corner);
  vtkQtChartAxis *yAxis = layer->getVerticalAxis(corner);

  int domainIndex = -1;
  const vtkQtChartSeriesDomain *seriesDomain =
      this->Internal->Domain.getDomain(xAxis->getAxisDomain(),
      yAxis->getAxisDomain(), &domainIndex);
  if(seriesDomain)
    {
    // Set up the painter clipping and offset for panning.
    painter->setClipRect(this->Internal->Bounds);
    painter->translate(-space->getXOffset(), -space->getYOffset());

    // Get the list of series in the selected domain.
    vtkQtStackedChartSeries *series = 0;
    vtkQtStackedChartSeriesOptions *options = 0;
    QList<int> seriesList = this->Internal->Groups.getGroup(domainIndex);
    QMutableListIterator<int> iter(seriesList);
    iter.toBack();
    while(iter.hasPrevious())
      {
      // Set up the painter for the series.
      int index = iter.previous();
      series = this->Internal->Series[index];
      options = this->getStackedSeriesOptions(index);
      QColor light = vtkQtChartColors::lighter(options->getBrush().color());
      painter->setPen(options->getPen());
      if(series->IsHighlighted)
        {
        painter->setBrush(light);
        }
      else if(this->Options->isGradientDislpayed())
        {
        QLinearGradient gradient(series->Gradient1, series->Gradient2);
        QColor color = options->getBrush().color();
        gradient.setColorAt(0.0, color);
        gradient.setColorAt(1.0, color.dark());
        painter->setBrush(QBrush(gradient));
        }
      else
        {
        painter->setBrush(options->getBrush());
        }

      // Draw the series polygon.
      painter->drawPolygon(*series->Polygon);

      // Draw the point highlights for the series.
      painter->setBrush(light);
      QList<QPolygonF *>::Iterator highlight = series->Highlights.begin();
      for( ; highlight != series->Highlights.end(); ++highlight)
        {
        painter->drawPolygon(*(*highlight));
        }
      }
    }
}

void vtkQtStackedChart::reset()
{
  // Make sure the selection model is notified of the change.
  this->InModelChange = true;
  this->Selection->beginModelReset();

  // Clean up the old view items.
  bool needsLayout = this->Internal->Series.size() > 0;
  QList<vtkQtStackedChartSeries *>::Iterator iter =
      this->Internal->Series.begin();
  for( ; iter != this->Internal->Series.end(); ++iter)
    {
    delete *iter;
    }

  this->Internal->Series.clear();
  this->Internal->Domain.clear();
  this->Internal->Groups.clear();

  // Add items for the new model.
  if(this->Model && this->ChartArea)
    {
    int total = this->Model->getNumberOfSeries();
    if(total > 0)
      {
      if(needsLayout)
        {
        needsLayout = false;
        emit this->rangeChanged();
        }

      this->insertSeries(0, total - 1);
      }
    }

  if(needsLayout)
    {
    emit this->rangeChanged();
    emit this->layoutNeeded();
    }

  // Notify the slection model that the reset is complete, which may
  // generate a selection changed signal.
  this->Selection->endModelReset();
  this->InModelChange = false;
}

vtkQtChartSeriesOptions *vtkQtStackedChart::createOptions(
    QObject *parentObject)
{
  return new vtkQtStackedChartSeriesOptions(parentObject);
}

void vtkQtStackedChart::setupOptions(vtkQtChartSeriesOptions *options)
{
  vtkQtStackedChartSeriesOptions *seriesOptions =
      qobject_cast<vtkQtStackedChartSeriesOptions *>(options);
  if(seriesOptions)
    {
    // Listen for series options changes.
    this->connect(seriesOptions, SIGNAL(visibilityChanged(bool)),
        this, SLOT(handleSeriesVisibilityChange(bool)));
    this->connect(seriesOptions, SIGNAL(penChanged(const QPen &)),
        this, SLOT(handleSeriesPenChange(const QPen &)));
    this->connect(seriesOptions, SIGNAL(brushChanged(const QBrush &)),
        this, SLOT(handleSeriesBrushChange(const QBrush &)));
    }
}

void vtkQtStackedChart::prepareSeriesInsert(int first, int last)
{
  if(this->ChartArea)
    {
    // Notify the selection model of the change. The selection will be
    // adjusted for the changes in this call so it can be layed out
    // when the changes are completed.
    this->InModelChange = true;
    this->Selection->beginInsertSeries(first, last);
    }
}

void vtkQtStackedChart::insertSeries(int first, int last)
{
  if(this->ChartArea)
    {
    // Update the series indexes stored in the domain groups.
    this->Internal->Groups.prepareInsert(first, last);

    // Add an item for each series.
    QList<int> tableGroups;
    vtkQtStackedChartSeriesOptions *options = 0;
    for(int i = first; i <= last; i++)
      {
      // Only add a polygon if the series y-axis range is numeric.
      QPolygonF *polygon = 0;
      QList<QVariant> yDomain = this->Model->getSeriesRange(i, 1);
      if(yDomain.size() == 2)
        {
        QVariant::Type domain = yDomain[0].type();
        if(domain == QVariant::Int || domain == QVariant::Double)
          {
          polygon = new QPolygonF();
          }
        }

      this->Internal->Series.insert(i, new vtkQtStackedChartSeries(polygon));
      options = this->getStackedSeriesOptions(i);
      if(polygon && options->isVisible())
        {
        // Add the series to the domain if it is visible.
        int seriesGroup = -1;
        this->addSeriesDomain(i, &seriesGroup);
        if(seriesGroup != -1 && !tableGroups.contains(seriesGroup))
          {
          tableGroups.append(seriesGroup);
          }
        }
      }

    if(tableGroups.size() > 0)
      {
      QList<int>::Iterator iter = tableGroups.begin();
      for( ; iter != tableGroups.end(); ++iter)
        {
        this->updateItemMap(*iter);
        this->createTable(*iter);
        this->createQuadTable(*iter);
        }

      emit this->rangeChanged();
      emit this->layoutNeeded();
      }

    // Close the event for the selection model, which will trigger a
    // selection change signal.
    this->Selection->endInsertSeries(first, last);
    this->InModelChange = false;
    }
}

void vtkQtStackedChart::startSeriesRemoval(int first, int last)
{
  if(this->ChartArea)
    {
    // Notify the selection model of the change. The selection will be
    // adjusted for the changes in this call so it can be layed out
    // when the changes are completed.
    this->InModelChange = true;
    this->Selection->beginRemoveSeries(first, last);

    // Remove each of the series items.
    for( ; last >= first; last--)
      {
      delete this->Internal->Series.takeAt(last);
      }
    }
}

void vtkQtStackedChart::finishSeriesRemoval(int first, int last)
{
  if(this->ChartArea)
    {
    // Find which groups need to be re-calculated
    QList<int> groups;
    QList<int>::Iterator iter;
    for(int i = first; i <= last; i++)
      {
      int index = this->Internal->Groups.removeSeries(i);
      if(index != -1)
        {
        // Add the group indexes in reverse order.
        bool doAdd = true;
        for(iter = groups.begin(); iter != groups.end(); ++iter)
          {
          if(index > *iter)
            {
            doAdd = false;
            groups.insert(iter, index);
            break;
            }
          else if(index == *iter)
            {
            doAdd = false;
            break;
            }
          }

        if(doAdd)
          {
          groups.append(index);
          }
        }
      }

    for(iter = groups.begin(); iter != groups.end(); ++iter)
      {
      if(this->Internal->Groups.getNumberOfSeries(*iter) == 0)
        {
        // Remove the empty domain.
        this->Internal->Domain.removeDomain(*iter);
        }
      else
        {
        // Re-calculate the chart domain and table.
        this->updateItemMap(*iter);
        this->calculateXDomain(*iter);
        this->createTable(*iter);
        this->createQuadTable(*iter);
        }
      }

    // Fix the stored indexes in the domain groups.
    this->Internal->Groups.finishRemoval(first, last);
    if(groups.size() > 0)
      {
      emit this->rangeChanged();
      emit this->layoutNeeded();
      }

    // Close the event for the selection model, which will trigger a
    // selection change signal.
    this->Selection->endRemoveSeries(first, last);
    this->InModelChange = false;
    }
}

void vtkQtStackedChart::handleAxesCornerChange()
{
  if(this->Model && this->ChartArea)
    {
    emit this->rangeChanged();
    emit this->layoutNeeded();
    }
}

void vtkQtStackedChart::handleSumationChange()
{
  if(this->Model && this->ChartArea)
    {
    for(int i = 0; i < this->Internal->Groups.getNumberOfGroups(); i++)
      {
      if(this->Options->isSumNormalized())
        {
        this->normalizeTable(i);
        this->calculateYDomain(i);
        }
      else
        {
        this->createTable(i);
        }
      }

    if(this->Internal->Groups.getNumberOfGroups() > 0)
      {
      emit this->rangeChanged();
      emit this->layoutNeeded();
      }
    }
}

void vtkQtStackedChart::handleGradientChange()
{
  if(this->Model && this->ChartArea)
    {
    if(this->Options->isGradientDislpayed())
      {
      // Update the gradient points.
      QList<vtkQtStackedChartSeries *>::Iterator iter =
          this->Internal->Series.begin();
      for(int i = 0; iter != this->Internal->Series.end(); ++iter, ++i)
        {
        if((*iter)->Polygon)
          {
          (*iter)->updateGradient();
          }
        }
      }

    this->update();
    }
}

void vtkQtStackedChart::handleSeriesVisibilityChange(bool visible)
{
  // Get the series index from the options index.
  vtkQtStackedChartSeriesOptions *options =
      qobject_cast<vtkQtStackedChartSeriesOptions *>(this->sender());
  int series = this->getSeriesOptionsIndex(options);
  if(series >= 0 && series < this->Internal->Series.size() &&
      this->Internal->Series[series]->Polygon)
    {
    if(visible)
      {
      int seriesGroup = -1;
      this->addSeriesDomain(series, &seriesGroup);
      if(seriesGroup != -1)
        {
        this->updateItemMap(seriesGroup);
        this->createTable(seriesGroup);
        this->createQuadTable(seriesGroup);
        emit this->rangeChanged();
        emit this->layoutNeeded();
        }
      }
    else
      {
      this->Internal->Series[series]->setMapping(-1, -1);
      int seriesGroup = this->Internal->Groups.removeSeries(series);
      if(seriesGroup != -1)
        {
        // If the group is empty, remove the domain.
        if(this->Internal->Groups.getNumberOfSeries(seriesGroup) == 0)
          {
          this->Internal->Domain.removeDomain(seriesGroup);
          }
        else
          {
          // Re-calculate the domain.
          this->updateItemMap(seriesGroup);
          this->calculateXDomain(seriesGroup);
          this->createTable(seriesGroup);
          this->createQuadTable(seriesGroup);
          }

        this->Internal->Groups.finishRemoval();
        emit this->rangeChanged();
        emit this->layoutNeeded();
        }
      }
    }
}

void vtkQtStackedChart::handleSeriesPenChange(const QPen &)
{
  // Get the series index from the options.
  vtkQtStackedChartSeriesOptions *options =
      qobject_cast<vtkQtStackedChartSeriesOptions *>(this->sender());
  int series = this->getSeriesOptionsIndex(options);
  if(series >= 0 && series < this->Internal->Series.size())
    {
    this->update();
    }
}

void vtkQtStackedChart::handleSeriesBrushChange(const QBrush &)
{
  // Get the series index from the options.
  vtkQtStackedChartSeriesOptions *options =
      qobject_cast<vtkQtStackedChartSeriesOptions *>(this->sender());
  int series = this->getSeriesOptionsIndex(options);
  if(series >= 0 && series < this->Internal->Series.size())
    {
    this->update();
    }
}

void vtkQtStackedChart::updateHighlights()
{
  if(!this->InModelChange && this->ChartArea)
    {
    // Remove the current selection.
    QList<vtkQtStackedChartSeries *>::Iterator iter =
        this->Internal->Series.begin();
    for( ; iter != this->Internal->Series.end(); ++iter)
      {
      (*iter)->IsHighlighted = false;
      (*iter)->clearHighlights();
      }

    // Get the current selection from the selection model.
    if(!this->Selection->isSelectionEmpty())
      {
      const vtkQtChartSeriesSelection &current =
          this->Selection->getSelection();
      if(current.getType() == vtkQtChartSeriesSelection::SeriesSelection)
        {
        const vtkQtChartIndexRangeList &series = current.getSeries();
        vtkQtChartIndexRangeList::ConstIterator jter = series.begin();
        for( ; jter != series.end(); ++jter)
          {
          for(int i = jter->first; i <= jter->second; i++)
            {
            this->Internal->Series[i]->IsHighlighted = true;
            }
          }
        }
      else if(current.getType() == vtkQtChartSeriesSelection::PointSelection)
        {
        this->layoutHighlights();
        }
      }

    this->update();
    }
}

void vtkQtStackedChart::layoutHighlights()
{
  if(this->Internal->Series.size() > 0 && !this->Selection->isSelectionEmpty())
    {
    // Get the current selection from the selection model.
    const vtkQtChartSeriesSelection &current =
        this->Selection->getSelection();
    if(current.getType() == vtkQtChartSeriesSelection::PointSelection)
      {
      const QList<vtkQtChartSeriesSelectionItem> &points =
          current.getPoints();
      QList<vtkQtChartSeriesSelectionItem>::ConstIterator jter;
      for(jter = points.begin(); jter != points.end(); ++jter)
        {
        // Clear the current highlights.
        vtkQtStackedChartSeries *item = this->Internal->Series[jter->Series];
        item->clearHighlights();

        // Add lightened polygons for the selected points.
        int half = item->Polygon->size() / 2;
        vtkQtChartIndexRangeList::ConstIterator kter = jter->Points.begin();
        for( ; kter != jter->Points.end(); ++kter)
          {
          // Add the mid-point to the front if needed.
          QPolygonF *selectedPoints = new QPolygonF();
          if(kter->first != 0)
            {
            selectedPoints->append(this->Internal->getMidPoint(
                (*item->Polygon)[kter->first - 1],
                (*item->Polygon)[kter->first]));
            }

          // Add the selected points.
          int count = kter->second - kter->first + 1;
          *selectedPoints << item->Polygon->mid(kter->first, count);

          // Add a midpoint to the end if needed. Add one for the
          // beginning of the bottom half as well.
          int bSecond = item->Polygon->size() - 1 - kter->first;
          int bFirst = bSecond - count + 1;
          if(kter->second < half - 1)
            {
            selectedPoints->append(this->Internal->getMidPoint(
                (*item->Polygon)[kter->second],
                (*item->Polygon)[kter->second + 1]));
            selectedPoints->append(this->Internal->getMidPoint(
                (*item->Polygon)[bFirst - 1], (*item->Polygon)[bFirst]));
            }

          // Add the selected points for the bottom half.
          *selectedPoints << item->Polygon->mid(bFirst, count);

          // Add the final mid-point if needed.
          if(kter->first != 0)
            {
            selectedPoints->append(this->Internal->getMidPoint(
                (*item->Polygon)[bSecond], (*item->Polygon)[bSecond + 1]));
            }

          // Add the highlight polygon.
          item->Highlights.append(selectedPoints);
          }
        }
      }
    }
}

void vtkQtStackedChart::addSeriesDomain(int series, int *seriesGroup)
{
  QList<QVariant> xDomain;
  QList<QVariant> yDomain = this->Model->getSeriesRange(series, 1);
  int points = this->Model->getNumberOfSeriesValues(series);
  for(int j = 0; j < points; j++)
    {
    xDomain.append(this->Model->getSeriesValue(series, j, 0));
    }

  // The y-axis domain is needed to separate the series groups.
  vtkQtChartSeriesDomain seriesDomain;
  seriesDomain.getXDomain().setDomain(xDomain);
  seriesDomain.getYDomain().setRange(yDomain);
  this->Internal->Domain.mergeDomain(seriesDomain, seriesGroup);

  // Add the series index to the domain group.
  this->Internal->Groups.insertSeries(series, *seriesGroup);
}

void vtkQtStackedChart::updateItemMap(int seriesGroup)
{
  QList<int> groupSeries = this->Internal->Groups.getGroup(seriesGroup);
  QList<int>::Iterator iter = groupSeries.begin();
  for(int i = 0; iter != groupSeries.end(); ++iter, ++i)
    {
    this->Internal->Series[*iter]->setMapping(seriesGroup, i);
    }
}

void vtkQtStackedChart::createTable(int seriesGroup)
{
  // Clear the group table and the associated y-axis domain.
  vtkQtStackedChartSeriesGroup *group =
      this->Internal->Groups.Tables[seriesGroup];
  group->Data.clear();
  vtkQtChartSeriesDomain *seriesDomain =
      this->Internal->Domain.getDomain(seriesGroup);
  seriesDomain->getYDomain().clear();

  // Get the x-axis domain.
  bool isRange = false;
  QList<QVariant> xDomain = seriesDomain->getXDomain().getDomain(isRange);
  if(xDomain.size() > 0)
    {
    // Get the list of series for the group.
    QList<int> seriesList = this->Internal->Groups.getGroup(seriesGroup);
    QList<int>::Iterator iter = seriesList.begin();
    for(int i = 0; iter != seriesList.end(); ++iter, ++i)
      {
      // Add a table row for each of the series.
      int k = 0;
      QVariant xValue, yValue;
      group->Data.append(QVector<double>(xDomain.size(), 0.0));
      int points = this->Model->getNumberOfSeriesValues(*iter);
      for(int j = 0; j < points; j++, k++)
        {
        // Find the matching x-axis value in the domain.
        xValue = this->Model->getSeriesValue(*iter, j, 0);
        while(k < xDomain.size() && xValue != xDomain[k])
          {
          if(i > 0)
            {
            group->Data[i][k] = group->Data[i - 1][k];
            }

          k++;
          }

        if(k >= xDomain.size())
          {
          break;
          }

        // Get the y-axis value.
        yValue = this->Model->getSeriesValue(*iter, j, 1);
        group->Data[i][k] = yValue.toDouble();

        // Stack the series by adding the previous series value.
        if(i > 0)
          {
          group->Data[i][k] += group->Data[i - 1][k];
          }
        }

      // Fill in any remaining table columns.
      if(i > 0)
        {
        for( ; k < xDomain.size(); k++)
          {
          group->Data[i][k] = group->Data[i - 1][k];
          }
        }
      }

    // Normalize the table if the user requested it.
    if(this->Options->isSumNormalized())
      {
      this->normalizeTable(seriesGroup);
      }

    this->calculateYDomain(seriesGroup);
    }
}

void vtkQtStackedChart::normalizeTable(int seriesGroup)
{
  vtkQtStackedChartSeriesGroup *group =
      this->Internal->Groups.Tables[seriesGroup];
  if(group->Data.size() == 0)
    {
    return;
    }

  int last = group->Data.size() - 1;
  int count = group->Data[0].size();
  for(int j = 0; j < count; j++)
    {
    double total = group->Data[last][j];
    if(total > 0)
      {
      int i = 0;
      for( ; i < group->Data.size(); i++)
        {
        double fraction = group->Data[i][j] / total;
        group->Data[i][j] = 100.0 * fraction;
        }
      }
    }
}

void vtkQtStackedChart::calculateXDomain(int seriesGroup)
{
  vtkQtChartSeriesDomain *seriesDomain =
      this->Internal->Domain.getDomain(seriesGroup);
  seriesDomain->getXDomain().clear();

  // Get the list of series in the group and merge the domains.
  QList<int> seriesList = this->Internal->Groups.getGroup(seriesGroup);
  QList<int>::Iterator iter = seriesList.begin();
  for( ; iter != seriesList.end(); ++iter)
    {
    QList<QVariant> xDomain;
    int points = this->Model->getNumberOfSeriesValues(*iter);
    for(int j = 0; j < points; j++)
      {
      xDomain.append(this->Model->getSeriesValue(*iter, j, 0));
      }

    seriesDomain->getXDomain().mergeDomain(xDomain);
    }
}

void vtkQtStackedChart::calculateYDomain(int seriesGroup)
{
  vtkQtStackedChartSeriesGroup *group =
      this->Internal->Groups.Tables[seriesGroup];
  vtkQtChartSeriesDomain *seriesDomain =
      this->Internal->Domain.getDomain(seriesGroup);
  seriesDomain->getYDomain().clear();

  // Use the first and last rows of the table to determine the minimum
  // and maximum respectively.
  if(group->Data.size() > 0)
    {
    double minimum = 0;
    double maximum = 0;
    QVector<double>::Iterator iter = group->Data[0].begin();
    QVector<double>::Iterator rowEnd = group->Data[0].end();
    QVector<double>::Iterator jter = group->Data.last().begin();
    if(iter != rowEnd)
      {
      minimum = *iter;
      maximum = *jter;
      ++iter;
      ++jter;
      }

    for( ; iter != rowEnd; ++iter, ++jter)
      {
      if(*iter < minimum)
        {
        minimum = *iter;
        }

      if(*jter > maximum)
        {
        maximum = *jter;
        }
      }

    QList<QVariant> yDomain;
    yDomain.append(QVariant(minimum));
    yDomain.append(QVariant(maximum));
    seriesDomain->getYDomain().setRange(yDomain);
    }
}

void vtkQtStackedChart::createQuadTable(int seriesGroup)
{
  // Clear the current quad table.
  vtkQtStackedChartSeriesGroup *group =
      this->Internal->Groups.Tables[seriesGroup];
  group->Shapes.clear();

  // Build the quad table from the value table sizes.
  int numSeries = group->Data.size();
  if(numSeries > 0)
    {
    int points = (group->Data[0].size() - 1) * 2;
    if(points > 0)
      {
      // Get the series list for the domain group.
      QList<int> seriesList = this->Internal->Groups.getGroup(seriesGroup);

      // Create the list of quads for each series in the group.
      int j = 0;
      vtkQtStackedChartSeries *series = 0;
      QList<int>::Iterator iter = seriesList.begin();
      for( ; iter != seriesList.end(); ++iter)
        {
        series = this->Internal->Series[*iter];
        series->clearQuads();
        for(j = 0; j < points; j++)
          {
          // Get the series index for the quad. There are two quads per
          // interval for selecting and highlighting points.
          int index = j / 2;
          if(j % 2 > 0)
            {
            index++;
            }

          series->Quads.append(new vtkQtChartQuad(*iter, index));
          }
        }

      for(j = 0; j < points; j++)
        {
        // Add a list for the y-direction quads.
        group->Shapes.append(QList<vtkQtChartShape *>());
        for(int i = numSeries - 1; i >= 0; i--)
          {
          // Add a quad for each of the series at this x-interval.
          series = this->Internal->Series[seriesList[i]];
          group->Shapes.last().append(series->Quads[j]);
          }
        }
      }
    }
}

void vtkQtStackedChart::buildQuadTree(int seriesGroup)
{
  this->BuildNeeded = false;
  if(seriesGroup == this->Internal->CurrentGroup)
    {
    this->Internal->QuadTree.update();
    }
  else
    {
    this->Internal->CurrentGroup = seriesGroup;
    vtkQtStackedChartSeriesGroup *group =
        this->Internal->Groups.Tables[seriesGroup];
    this->Internal->QuadTree.build(group->Shapes);
    }
}


