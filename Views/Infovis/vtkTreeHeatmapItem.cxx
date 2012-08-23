/*=========================================================================

  Program:   Visualization Toolkit
  Module:    TestDiagram.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkTreeHeatmapItem.h"

#include "vtkBrush.h"
#include "vtkContext2D.h"
#include "vtkDataSetAttributes.h"
#include "vtkDoubleArray.h"
#include "vtkGraphLayout.h"
#include "vtkImageData.h"
#include "vtkLookupTable.h"
#include "vtkMarkerUtilities.h"
#include "vtkNew.h"
#include "vtkObjectFactory.h"
#include "vtkPen.h"
#include "vtkStringArray.h"
#include "vtkTable.h"
#include "vtkTextProperty.h"
#include "vtkTree.h"
#include "vtkTreeLayoutStrategy.h"

#include <float.h>
#include <sstream>

vtkStandardNewMacro(vtkTreeHeatmapItem);
vtkCxxSetObjectMacro(vtkTreeHeatmapItem, Tree, vtkTree);
vtkCxxSetObjectMacro(vtkTreeHeatmapItem, Table, vtkTable);

vtkTreeHeatmapItem::vtkTreeHeatmapItem()
{
  this->TreeHeatmapBuildTime = 0;
  this->Tree = 0;
  this->Table = 0;
  this->Multiplier = 100.0;
}

vtkTreeHeatmapItem::~vtkTreeHeatmapItem()
{
  if (this->Tree)
    {
    this->Tree->Delete();
    }

  if (this->Table)
    {
    this->Table->Delete();
    }

  while (!this->LookupTables.empty())
    {
    this->LookupTables.back()->Delete();
    this->LookupTables.pop_back();
    }
}

bool vtkTreeHeatmapItem::Paint(vtkContext2D *painter)
{
  if (!this->Tree || !this->Table)
    {
    return true;
    }

  if (this->IsDirty())
    {
    this->RebuildBuffers();
    }

  this->PaintBuffers(painter);
  return true;
}

bool vtkTreeHeatmapItem::IsDirty()
{
  if (!this->Tree || !this->Table)
    {
    return false;
    }
  if (this->Tree->GetMTime() > this->TreeHeatmapBuildTime ||
      this->Table->GetMTime() > this->GetMTime())
    {
    return true;
    }
  return false;
}

void vtkTreeHeatmapItem::RebuildBuffers()
{
  vtkNew<vtkTreeLayoutStrategy> strategy;
  strategy->SetDistanceArrayName("node weight");
  strategy->SetLeafSpacing(1.0);
  strategy->SetRotation(90.0);

  this->Layout->SetLayoutStrategy(strategy.GetPointer());
  this->Layout->SetInputData(this->Tree);
  this->Layout->Update();
  this->LayoutTree = vtkTree::SafeDownCast(this->Layout->GetOutput());

  this->ComputeMultiplier();
  this->InitializeLookupTables();

  if(this->Tree->GetMTime() > this->Table->GetMTime())
    {
    this->TreeHeatmapBuildTime = this->Tree->GetMTime();
    }
  else
    {
    this->TreeHeatmapBuildTime = this->Table->GetMTime();
    }
}
  
// Figure out the multiplier we need to use so that the table cells are large
// enough to be labeled.  Currently, we require the boxes to be big enough so
// an 18 pt font can be used.
void vtkTreeHeatmapItem::ComputeMultiplier()
{
  double targetFontSize = 18;
  double yMax = DBL_MIN; 
  double targetPoint[3];
  for (vtkIdType edge = 0; edge < this->LayoutTree->GetNumberOfEdges(); ++edge)
    {
    vtkIdType target = this->LayoutTree->GetTargetVertex(edge);
    this->LayoutTree->GetPoint(target, targetPoint);
    if (targetPoint[1] > yMax)
      {
      yMax = targetPoint[1];
      }
    }
  double currentFontSize =
    (yMax * this->Multiplier) / this->Table->GetNumberOfRows();
  if (currentFontSize < targetFontSize)
    {
    this->Multiplier = (this->Table->GetNumberOfRows() * targetFontSize) / yMax;
    }
}

void vtkTreeHeatmapItem::InitializeLookupTables()
{
  while (!this->LookupTables.empty())
    {
    this->LookupTables.back()->Delete();
    this->LookupTables.pop_back();
    }
  this->LookupTables.reserve( this->Table->GetNumberOfColumns() + 1 );

  for (vtkIdType column = 1; column < this->Table->GetNumberOfColumns();
       ++column)
    {
    if (this->Table->GetValue(0, column).IsString())
      {
      this->GenerateLookupTableForStringColumn(column);
      continue;
      }
    double min = DBL_MAX;
    double max = DBL_MIN;
    vtkLookupTable *lookupTable = vtkLookupTable::New();
    this->LookupTables.push_back(lookupTable);
    for (vtkIdType row = 0; row < this->Table->GetNumberOfRows(); ++row)
      {
      double value = this->Table->GetValue(row, column).ToDouble();
      if (value > max)
        {
        max = value;
        }
      else if (value < min)
        {
        min = value;
        }
      }
    this->LookupTables[column-1]->SetNumberOfTableValues(256);
    this->LookupTables[column-1]->SetRange(min, max);
    this->LookupTables[column-1]->Build();
    }
}
      
void vtkTreeHeatmapItem::GenerateLookupTableForStringColumn(vtkIdType column)
{
  //generate a sorted vector of all the strings in this column
  std::vector< std::string > sortedStrings;
  for (vtkIdType row = 0; row < this->Table->GetNumberOfRows(); ++row)
    {
    std::string value = this->Table->GetValue(row, column).ToString();
    // no duplicates
    if (std::find(sortedStrings.begin(), sortedStrings.end(), value) ==
        sortedStrings.end())
      {
      sortedStrings.push_back(value);
      }
    }
  std::sort(sortedStrings.begin(), sortedStrings.end());

  // map each string to a double value based on its position
  // in alphabetical order
  std::map< std::string, double> stringToDouble;
  for (int i = 0; i < sortedStrings.size(); ++i)
    {
    stringToDouble[ sortedStrings[i] ] = (double)i;
    }
 
  // associate this mapping with the column number
  this->StringToDoubleMaps[column] = stringToDouble;
   
  // generate a lookup table for this column
  this->LookupTables.push_back(vtkLookupTable::New());
  this->LookupTables[column-1]->SetNumberOfTableValues(256);
  this->LookupTables[column-1]->SetRange(0, sortedStrings.size() - 1);
  this->LookupTables[column-1]->Build();
}

void vtkTreeHeatmapItem::PaintBuffers(vtkContext2D *painter)
{
  double yMax = DBL_MIN; 
  double xMax = DBL_MIN; 
  double sourcePoint[3];
  double targetPoint[3];

  for (vtkIdType edge = 0; edge < this->LayoutTree->GetNumberOfEdges(); ++edge)
    {
    vtkIdType source = this->LayoutTree->GetSourceVertex(edge);
    vtkIdType target = this->LayoutTree->GetTargetVertex(edge);

    this->LayoutTree->GetPoint(source, sourcePoint);
    this->LayoutTree->GetPoint(target, targetPoint);

    double x0 = sourcePoint[0] * this->Multiplier;
    double y0 = sourcePoint[1] * this->Multiplier;
    double x1 = targetPoint[0] * this->Multiplier; 
    double y1 = targetPoint[1] * this->Multiplier;

    painter->DrawLine (x0, y0, x0, y1);
    painter->DrawLine (x0, y1, x1, y1);

    if (x1 > xMax)
      {
      xMax = x1;
      }
    if (y1 > yMax)
      {
      yMax = y1;
      }
    }

  // calculate how large our table cells will be when they are drawn
  double cellHeight = yMax / this->Table->GetNumberOfRows();
  double cellWidth; 
  if (cellHeight * 2 > 100)
    {
    cellWidth = 100;
    }
  else
    {
    cellWidth = cellHeight * 2;
    }
  
  // leave a small amount of space between the tree, the table,
  // and the row/column labels
  double spacing = cellWidth * 0.25;

  // get array of node names from the tree
  vtkStringArray *nodeNames = vtkStringArray::SafeDownCast(
    this->LayoutTree->GetVertexData()->GetAbstractArray("node name"));

  // get array of row names from the table.  We assume this is the first row.
  vtkStringArray *tableNames = vtkStringArray::SafeDownCast(
    this->Table->GetColumn(0));

  // set up our text property to draw row names
  painter->GetTextProp()->SetColor(0.0, 0.0, 0.0);
  painter->GetTextProp()->SetJustificationToLeft();
  painter->GetTextProp()->SetVerticalJustificationToCentered();
  painter->GetTextProp()->SetOrientation(0);
  
  // calculate a font size that's appropriate for this zoom level
  float stringBounds[4];
  stringBounds[3] = FLT_MAX;
  std::string testString = "Igq"; //selected for range of height
  int currentFontSize = floor(cellHeight);
  painter->GetTextProp()->SetFontSize(currentFontSize);
  painter->ComputeStringBounds(testString, stringBounds);
  if (stringBounds[3] > cellHeight)
    {
    while (stringBounds[3] > cellHeight && currentFontSize > 0)
      {
      --currentFontSize;
      painter->GetTextProp()->SetFontSize(currentFontSize);
      painter->ComputeStringBounds(testString, stringBounds);
      }
    }
  else
    {
      while (stringBounds[3] < cellHeight)
      {
      ++currentFontSize;
      painter->GetTextProp()->SetFontSize(currentFontSize);
      painter->ComputeStringBounds(testString, stringBounds);
      }
    --currentFontSize;
    painter->GetTextProp()->SetFontSize(currentFontSize);
    }
  
  double xStart, yStart;
  double yMaxTable = DBL_MIN;
  
  for (vtkIdType vertex = 0; vertex < this->LayoutTree->GetNumberOfVertices();
       ++vertex)
    {
    if (!this->LayoutTree->IsLeaf(vertex))
      {
      continue;
      }


    // find the row in the table that corresponds to this vertex
    double point[3];
    this->LayoutTree->GetPoint(vertex, point);
    std::string nodeName = nodeNames->GetValue(vertex);
    vtkIdType tableRow = tableNames->LookupValue(nodeName);

    for (vtkIdType column = 1; column < this->Table->GetNumberOfColumns();
         ++column)
      {
      // get the color for this cell from the lookup table
      double color[3];
      if (this->Table->GetValue(tableRow, column).IsString())
        {
        // get the string to int mapping for this column
        std::map< std::string, double> stringToDouble = this->StringToDoubleMaps[column];

        // get the integer value for the current string
        std::string cellStr = this->Table->GetValue(tableRow, column).ToString();
        double colorKey = stringToDouble[cellStr];
        
        // now we can lookup the appropriate color for this string
        this->LookupTables[column-1]->GetColor(colorKey, color);
        }
      else
        {
        vtkVariant value = this->Table->GetValue(tableRow, column);
        this->LookupTables[column-1]->GetColor(value.ToDouble(), color);
        }
      painter->GetBrush()->SetColorF(color[0], color[1], color[2]);

      // draw this cell of the table
      xStart = xMax + spacing + cellWidth * (column - 1);
      yStart = point[1] * this->Multiplier - (cellHeight / 2);
      painter->DrawRect(xStart, yStart, cellWidth, cellHeight);

      // keep track of where the top of the table is, so we know where to 
      // draw the column labels later.
      if (yStart + cellHeight > yMaxTable)
        {
        yMaxTable = yStart + cellHeight;
        }
      }

    // draw the label for this row
    xStart = xMax + spacing * 2 +
      cellWidth * (this->Table->GetNumberOfColumns() - 1);
    yStart = point[1] * this->Multiplier;
    painter->DrawString(xStart, yStart, nodeName);
    }
  
  // draw column labels
  painter->GetTextProp()->SetOrientation(90);
  for (vtkIdType column = 1; column < this->Table->GetNumberOfColumns();
       ++column)
    {
      std::string columnName = this->Table->GetColumn(column)->GetName();
      xStart = xMax + spacing + cellWidth * column - cellWidth / 2;
      yStart = yMaxTable + spacing;
      painter->DrawString(xStart, yStart, columnName);
    }
}

void vtkTreeHeatmapItem::PrintSelf(ostream &os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
  os << "Tree: " << (this->Tree ? "" : "(null)") << std::endl;
  if (this->Tree)
    {
    this->Tree->PrintSelf(os, indent.GetNextIndent());
    }
  os << "Table: " << (this->Table ? "" : "(null)") << std::endl;
  if (this->Table)
    {
    this->Table->PrintSelf(os, indent.GetNextIndent());
    }
}
