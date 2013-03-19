/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkGeoJSONWriter.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkGeoJSONWriter.h"

#include "vtkCellArray.h"
#include "vtkInformation.h"
#include "vtkLookupTable.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"

#include <vtksys/ios/sstream>

vtkStandardNewMacro(vtkGeoJSONWriter);

//------------------------------------------------------------------------------
vtkGeoJSONWriter::vtkGeoJSONWriter()
{
  this->FileName = NULL;
  this->OutputString = NULL;
  this->SetNumberOfOutputPorts(0);
  this->WriteToOutputString = false;
}

//------------------------------------------------------------------------------
vtkGeoJSONWriter::~vtkGeoJSONWriter()
{
  this->SetFileName(NULL);
  delete[] this->OutputString;
}

//------------------------------------------------------------------------------
void vtkGeoJSONWriter::PrintSelf(ostream & os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
  os << indent << "FileName: "
     << (this->FileName?this->FileName:"NONE") << endl;
  os << indent << "WriteToOutputString: "
     << (this->WriteToOutputString?"True":"False") << endl;
}

//------------------------------------------------------------------------------
int vtkGeoJSONWriter::FillInputPortInformation(int port, vtkInformation *info)
{
  if (port == 0)
    {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
    }
  return 1;
}

//------------------------------------------------------------------------------
ostream *vtkGeoJSONWriter::OpenFile()
{
  vtkDebugMacro(<<"Opening file\n");

  ostream *fptr;

  if (!this->WriteToOutputString)
    {
    if (!this->FileName)
      {
      vtkErrorMacro(<< "No FileName specified! Can't write!");
      return NULL;
      }

    fptr = new ofstream(this->FileName, ios::out);
    }
  else
    {
    // Get rid of any old output string.
    if (this->OutputString)
      {
      delete [] this->OutputString;
      this->OutputString = NULL;
      this->OutputStringLength = 0;
      }
    fptr = new vtksys_ios::ostringstream;
    }

  if (fptr->fail())
    {
    vtkErrorMacro(<< "Unable to open file: "<< this->FileName);
    delete fptr;
    return NULL;
    }

  return fptr;
}

//------------------------------------------------------------------------------
void vtkGeoJSONWriter::CloseFile(ostream *fp)
{
  vtkDebugMacro(<<"Closing file\n");

  if ( fp != NULL )
    {
    if (this->WriteToOutputString)
      {
      vtksys_ios::ostringstream *ostr =
        static_cast<vtksys_ios::ostringstream*>(fp);

      delete [] this->OutputString;
      this->OutputStringLength = static_cast<int>(ostr->str().size());
      //+1's account for null terminator
      this->OutputString = new char[ostr->str().size()+1];
      memcpy(this->OutputString, ostr->str().c_str(),
        this->OutputStringLength+1);
      }

    delete fp;
    }
}

//------------------------------------------------------------------------------
void vtkGeoJSONWriter::ConditionalComma(ostream *fp,
  vtkIdType cnt, vtkIdType limit)
{
  if (cnt+1 != limit)
    {
    *fp << ",";
    }
}

//------------------------------------------------------------------------------
void vtkGeoJSONWriter::WritePointColor(ostream *fp,
  vtkDataArray *da, vtkIdType ptId)
{
  if (da)
    {
    vtkLookupTable *lut = da->GetLookupTable();
    double b = da->GetTuple1(ptId);
    unsigned char *color = lut->MapValue(b);
    *fp << ","
        << (double)color[0]/255.0 << ","
        << (double)color[1]/255.0 << ","
        << (double)color[2]/255;
    }
}

//------------------------------------------------------------------------------
void vtkGeoJSONWriter::WriteData()
{
  ostream *fp;
  vtkPolyData *input = vtkPolyData::SafeDownCast(this->GetInput());

  vtkDebugMacro(<<"Writing vtk polygonal data to geojson file...");
  fp=this->OpenFile();
  if ( !fp )
    {
    return;
    }

  *fp << "{\n";
  *fp << "\"type\": \"Feature\",\n";
  *fp << "\"properties\": null,\n";

  *fp << "\"geometry\":\n";
  *fp << "{\n";
  *fp << "\"type\": \"GeometryCollection\",\n";
  *fp << "\"geometries\":\n";
  *fp << "[\n";

  //
  // Write polygonal data specific stuff
  //
  vtkIdType cellLoc = 0;
  vtkIdType *cellPts = NULL;
  vtkIdType cellSize = 0;
  vtkIdType numlines, numpolys;
  numlines = input->GetLines()->GetNumberOfCells();
  numpolys = input->GetPolys()->GetNumberOfCells();

  vtkDataArray *da = input->GetPointData()->GetScalars();
  if (!da)
    {
    da = input->GetPointData()->GetArray(0);
    }
  if (da)
    {
    vtkLookupTable *lut = da->GetLookupTable();
    if (!lut)
      {
      da->CreateDefaultLookupTable();
      lut = da->GetLookupTable();
      }
    lut->SetRange(da->GetRange());
    }
  for (int i = 0; i < 3; i++)
    {
    vtkCellArray *ca;
    switch (i)
      {
      case 0:
        ca = input->GetVerts();
        break;
      case 1:
        ca = input->GetLines();
        break;
      default:
        ca = input->GetPolys();
      }
    if (ca && ca->GetNumberOfCells())
      {
      switch (i)
        {
        case 0:
          *fp << "{\n";
          *fp << "\"type\": \"MultiPoint\",\n";
          break;
        case 1:
          *fp << "{\n";
          *fp << "\"type\": \"MultiLineString\",\n";
          break;
        default:
          *fp << "{\n";
          *fp << "\"type\": \"MultiPolygon\",\n";
        }
      *fp << "\"coordinates\":\n";
      *fp << "[\n";
      vtkIdType inCell;

      switch (i)
        {
        case 0:
          for (inCell = 0; inCell < ca->GetNumberOfCells(); inCell++)
            {
            ca->GetCell(cellLoc, cellSize, cellPts);
            cellLoc += cellSize+1;
            vtkIdType inPt;
            for (inPt = 0; inPt < cellSize; inPt++)
              {
              double coords[3];
              input->GetPoint(cellPts[inPt], coords);
              *fp << "[" << coords[0] << "," << coords[1] << "," << coords[2];
              this->WritePointColor(fp, da, cellPts[inPt]);
              *fp << "]";
              this->ConditionalComma(fp, inPt, cellSize);
              }
            this->ConditionalComma(fp, inCell, ca->GetNumberOfCells());
            *fp << "\n";
            }
          break;
        case 1:
          for (inCell = 0; inCell < ca->GetNumberOfCells(); inCell++)
            {
            *fp << "[ "; //one cell
            ca->GetCell(cellLoc, cellSize, cellPts);
            cellLoc += cellSize+1;
            vtkIdType inPt;
            for (inPt = 0; inPt < cellSize; inPt++)
              {
              double coords[3];
              input->GetPoint(cellPts[inPt], coords);
              *fp << "[" << coords[0] << "," << coords[1] << "," << coords[2];
              this->WritePointColor(fp, da, cellPts[inPt]);
              *fp << "]";
              this->ConditionalComma(fp, inPt, cellSize);
              }
           *fp << " ]";//one cell
           this->ConditionalComma(fp, inCell, ca->GetNumberOfCells());
           *fp << "\n";
           }
           break;
        default:
          for (inCell = 0; inCell < ca->GetNumberOfCells(); inCell++)
            {
            *fp << "[[ "; //one cell
            ca->GetCell(cellLoc, cellSize, cellPts);
            cellLoc += cellSize+1;
            vtkIdType inPt;
            for (inPt = 0; inPt < cellSize; inPt++)
              {
              double coords[3];
              input->GetPoint(cellPts[inPt], coords);
              *fp << "[" << coords[0] << "," << coords[1] << "," << coords[2];
              this->WritePointColor(fp, da, cellPts[inPt]);
              *fp << "]";
              this->ConditionalComma(fp, inPt, cellSize);
              }
            *fp << " ]]";//one cell
            this->ConditionalComma(fp, inCell, ca->GetNumberOfCells());
            *fp << "\n";
            }
        }
      *fp << "]"; //coordinates for this cell array
      if ((i == 0 && (numlines || numpolys)) ||
        (i == 1 && numpolys))
        {
        *fp << ",";
        }
      *fp << "\n";
      *fp << "}\n"; //this cell array
      }
    }
  *fp << "]\n";//feature.geometry.GeometryCollection.geometries
  *fp << "}\n";//feature.geometry
  *fp << "}\n";//feature

  fp->flush();
  if (fp->fail())
    {
    vtkErrorMacro("Problem writing result check disk space.");
    delete fp;
    fp = NULL;
    }

  this->CloseFile(fp);
}

//------------------------------------------------------------------------------
char *vtkGeoJSONWriter::RegisterAndGetOutputString()
{
  char *tmp = this->OutputString;

  this->OutputString = NULL;
  this->OutputStringLength = 0;

  return tmp;
}

//------------------------------------------------------------------------------
vtkStdString vtkGeoJSONWriter::GetOutputStdString()
{
  return vtkStdString(this->OutputString, this->OutputStringLength);
}
