/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkArrayReader.cxx

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

#include "vtkArrayReader.h"

#include <vtkCommand.h>
#include <vtkDenseArray.h>
#include <vtkObjectFactory.h>
#include <vtkSmartPointer.h>
#include <vtkSparseArray.h>
#include <vtkUnicodeString.h>

#include <vtksys/ios/sstream>
#include <vtkstd/stdexcept>
#include <vtkstd/string>

vtkStandardNewMacro(vtkArrayReader);

namespace {

template<typename ValueT>
void ExtractValue(istream& stream, ValueT& value)
{
  stream >> value;
}

void ExtractValue(istream& stream, vtkStdString& value)
{
  vtkstd::getline(stream, value);
  vtkStdString::size_type begin, end;
  begin = 0;  end = value.size();
  while ((begin < end) && isspace(value[begin])) begin++;
  while ((begin < end) && isspace(value[end-1])) end--;
  value = value.substr(begin, end);
}

void ExtractValue(istream& stream, vtkUnicodeString& value)
{
  vtkstd::string buffer;
  vtkstd::getline(stream, buffer);
  vtkStdString::size_type begin, end;
  begin = 0;  end = buffer.size();
  while ((begin < end) && isspace(buffer[begin])) begin++;
  while ((begin < end) && isspace(buffer[end-1])) end--;
  buffer = buffer.substr(begin, end);
  value = vtkUnicodeString::from_utf8(buffer);
}

void ReadHeader(istream& stream, vtkArrayExtents& extents, vtkArrayExtents::SizeT& non_null_size, vtkArray* array)
{
  if(!array)
    throw vtkstd::runtime_error("Missing array.");

  // Load the array name ...
  vtkstd::string name;
  vtkstd::getline(stream, name);
  array->SetName(name);

  // Load array extents ...
  vtkstd::string extents_string;
  vtkstd::getline(stream, extents_string);
  vtkstd::istringstream extents_buffer(extents_string);

  vtkArrayExtents::CoordinateT extent;
  vtkstd::vector<vtkArrayExtents::CoordinateT> temp_extents;
  for(extents_buffer >> extent; extents_buffer; extents_buffer >> extent)
    temp_extents.push_back(extent);

  extents.SetDimensions(0);
  while(temp_extents.size() > 1)
    {
    const vtkArrayExtents::CoordinateT begin = temp_extents.front();
    temp_extents.erase(temp_extents.begin());
    const vtkArrayExtents::CoordinateT end = temp_extents.front();
    temp_extents.erase(temp_extents.begin());
    extents.Append(vtkArrayRange(begin, end));
    }

  if(extents.GetDimensions() < 1)
    throw vtkstd::runtime_error("Array cannot have fewer than one dimension.");

  if(temp_extents.empty())
    throw vtkstd::runtime_error("Missing non null size.");

  non_null_size = temp_extents.back();

  array->Resize(extents);

  // Load dimension-labels ...
  for(vtkArrayExtents::DimensionT i = 0; i != extents.GetDimensions(); ++i)
    {
    vtkstd::string label;
    vtkstd::getline(stream, label);
    array->SetDimensionLabel(i, label);
    }
}

void ReadEndianOrderMark(istream& stream, bool& swap_endian)
{
  // Load the endian-order mark ...
  vtkTypeUInt32 endian_order = 0;
  stream.read(
    reinterpret_cast<char*>(&endian_order),
    sizeof(endian_order));

  swap_endian = endian_order == 0x12345678 ? false : true;
}

template<typename ValueT>
vtkSparseArray<ValueT>* ReadSparseArrayBinary(istream& stream)
{
  // Create the array ...
  vtkSmartPointer<vtkSparseArray<ValueT> > array = vtkSmartPointer<vtkSparseArray<ValueT> >::New();

  // Read the file header ...
  vtkArrayExtents extents;
  vtkArrayExtents::SizeT non_null_size = 0;
  bool swap_endian = false;
  ReadHeader(stream, extents, non_null_size, array);
  ReadEndianOrderMark(stream, swap_endian);

  // Read the array NULL value ...
  ValueT null_value;
  stream.read(
    reinterpret_cast<char*>(&null_value),
    sizeof(ValueT));
  array->SetNullValue(null_value);

  // Read array coordinates ...
  array->ReserveStorage(non_null_size);

  for(vtkArray::DimensionT i = 0; i != array->GetDimensions(); ++i)
    {
    stream.read(
      reinterpret_cast<char*>(array->GetCoordinateStorage(i)),
      non_null_size * sizeof(vtkArray::CoordinateT));
    }

  // Read array values ...
  stream.read(
    reinterpret_cast<char*>(array->GetValueStorage()),
    non_null_size * sizeof(ValueT));

  array->Register(0);
  return array;
}

template<>
vtkSparseArray<vtkStdString>* ReadSparseArrayBinary<vtkStdString>(istream& stream)
{
  // Create the array ...
  vtkSmartPointer<vtkSparseArray<vtkStdString> > array = vtkSmartPointer<vtkSparseArray<vtkStdString> >::New();

  // Read the file header ...
  vtkArrayExtents extents;
  vtkArrayExtents::SizeT non_null_size = 0;
  bool swap_endian = false;
  ReadHeader(stream, extents, non_null_size, array);
  ReadEndianOrderMark(stream, swap_endian);

  // Read the array NULL value ...
  vtkstd::string null_value;
  for(int character = stream.get(); stream; character = stream.get())
    {
    if(character == 0)
      {
      array->SetNullValue(null_value);
      break;
      }
    else
      {
      null_value += character;
      }
    }

  // Read array coordinates ...
  array->ReserveStorage(non_null_size);

  for(vtkArray::DimensionT i = 0; i != array->GetDimensions(); ++i)
    {
    stream.read(
      reinterpret_cast<char*>(array->GetCoordinateStorage(i)),
      non_null_size * sizeof(vtkArray::CoordinateT));
    }

  // Read array values ...
  vtkstd::string buffer;
  vtkArray::SizeT n = 0;
  for(int character = stream.get(); stream; character = stream.get())
    {
    if(character == 0)
      {
      array->SetValueN(n++, buffer);
      buffer.resize(0);
      }
    else
      {
      buffer += character;
      }
    }

  array->Register(0);
  return array;
}

template<>
vtkSparseArray<vtkUnicodeString>* ReadSparseArrayBinary<vtkUnicodeString>(istream& stream)
{
  // Create the array ...
  vtkSmartPointer<vtkSparseArray<vtkUnicodeString> > array = vtkSmartPointer<vtkSparseArray<vtkUnicodeString> >::New();

  // Read the file header ...
  vtkArrayExtents extents;
  vtkArrayExtents::SizeT non_null_size = 0;
  bool swap_endian = false;
  ReadHeader(stream, extents, non_null_size, array);
  ReadEndianOrderMark(stream, swap_endian);

  // Read the array NULL value ...
  vtkstd::string null_value;
  for(int character = stream.get(); stream; character = stream.get())
    {
    if(character == 0)
      {
      array->SetNullValue(vtkUnicodeString::from_utf8(null_value));
      break;
      }
    else
      {
      null_value += character;
      }
    }

  // Read array coordinates ...
  array->ReserveStorage(non_null_size);

  for(vtkArray::DimensionT i = 0; i != array->GetDimensions(); ++i)
    {
    stream.read(
      reinterpret_cast<char*>(array->GetCoordinateStorage(i)),
      non_null_size * sizeof(vtkArray::CoordinateT));
    }

  // Read array values ...
  vtkstd::string buffer;
  vtkArray::SizeT n = 0;
  for(int character = stream.get(); stream; character = stream.get())
    {
    if(character == 0)
      {
      array->SetValueN(n++, vtkUnicodeString::from_utf8(buffer));
      buffer.resize(0);
      }
    else
      {
      buffer += character;
      }
    }

  array->Register(0);
  return array;
}

template<typename ValueT>
vtkDenseArray<ValueT>* ReadDenseArrayBinary(istream& stream)
{
  // Create the array ...
  vtkSmartPointer<vtkDenseArray<ValueT> > array = vtkSmartPointer<vtkDenseArray<ValueT> >::New();

  // Read the file header ...
  vtkArrayExtents extents;
  vtkArrayExtents::SizeT non_null_size = 0;
  bool swap_endian = false;
  ReadHeader(stream, extents, non_null_size, array);
  ReadEndianOrderMark(stream, swap_endian);

  // Read array values ...
  stream.read(
    reinterpret_cast<char*>(array->GetStorage()),
    non_null_size * sizeof(ValueT));

  array->Register(0);
  return array;
}

template<>
vtkDenseArray<vtkStdString>* ReadDenseArrayBinary<vtkStdString>(istream& stream)
{
  // Create the array ...
  vtkSmartPointer<vtkDenseArray<vtkStdString> > array = vtkSmartPointer<vtkDenseArray<vtkStdString> >::New();

  // Read the file header ...
  vtkArrayExtents extents;
  vtkArrayExtents::SizeT non_null_size = 0;
  bool swap_endian = false;
  ReadHeader(stream, extents, non_null_size, array);
  ReadEndianOrderMark(stream, swap_endian);

  // Read array values ...
  vtkstd::string buffer;
  vtkArray::SizeT n = 0;
  for(int character = stream.get(); stream; character = stream.get())
    {
    if(character == 0)
      {
      array->SetValueN(n++, buffer);
      buffer.resize(0);
      }
    else
      {
      buffer += character;
      }
    }

  array->Register(0);
  return array;
}

template<>
vtkDenseArray<vtkUnicodeString>* ReadDenseArrayBinary<vtkUnicodeString>(istream& stream)
{
  // Create the array ...
  vtkSmartPointer<vtkDenseArray<vtkUnicodeString> > array = vtkSmartPointer<vtkDenseArray<vtkUnicodeString> >::New();

  // Read the file header ...
  vtkArrayExtents extents;
  vtkArrayExtents::SizeT non_null_size = 0;
  bool swap_endian = false;
  ReadHeader(stream, extents, non_null_size, array);
  ReadEndianOrderMark(stream, swap_endian);

  // Read array values ...
  vtkstd::string buffer;
  vtkArray::SizeT n = 0;
  for(int character = stream.get(); stream; character = stream.get())
    {
    if(character == 0)
      {
      array->SetValueN(n++, vtkUnicodeString::from_utf8(buffer));
      buffer.resize(0);
      }
    else
      {
      buffer += character;
      }
    }

  array->Register(0);
  return array;
}

template<typename ValueT>
vtkSparseArray<ValueT>* ReadSparseArrayAscii(istream& stream)
{
  // Create the array ...
  vtkSmartPointer<vtkSparseArray<ValueT> > array = vtkSmartPointer<vtkSparseArray<ValueT> >::New();

  // Read the stream header ...
  vtkArrayExtents extents;
  vtkArrayExtents::SizeT non_null_size = 0;
  ReadHeader(stream, extents, non_null_size, array);

  if(non_null_size > extents.GetSize())
    throw vtkstd::runtime_error("Too many values for a sparse array.");

  // Read the array NULL value ...
  vtkstd::string line_buffer;
  vtkstd::getline(stream, line_buffer);
  if(!stream)
    throw vtkstd::runtime_error("Premature end-of-stream reading NULL value.");

  vtkstd::istringstream line_stream(line_buffer);
  ValueT null_value;
  ExtractValue(line_stream, null_value);
  if(!line_stream)
    throw vtkstd::runtime_error("Missing NULL value.");
  array->SetNullValue(null_value);

  // Setup storage for the stream contents ...
  array->ReserveStorage(non_null_size);
  vtkstd::vector<vtkArray::CoordinateT*> coordinates(array->GetDimensions());
  for(vtkArray::DimensionT j = 0; j != array->GetDimensions(); ++j)
    coordinates[j] = array->GetCoordinateStorage(j);
  ValueT* value = array->GetValueStorage();

  // Read the stream contents ...
  vtkArray::SizeT value_count = 0;
  for(vtkstd::getline(stream, line_buffer); stream; vtkstd::getline(stream, line_buffer), ++value_count)
    {
    if(value_count + 1 > non_null_size)
      throw vtkstd::runtime_error("Stream contains too many values.");

    line_stream.clear();
    line_stream.str(line_buffer);

    for(vtkArray::DimensionT j = 0; j != array->GetDimensions(); ++j)
      {
      line_stream >> *(coordinates[j] + value_count);
      if(!extents[j].Contains(*(coordinates[j] + value_count)))
        throw vtkstd::runtime_error("Coordinate out-of-bounds.");
      if(!line_stream)
        throw vtkstd::runtime_error("Missing coordinate.");
      }

    ExtractValue(line_stream, *(value + value_count));
    if(!line_stream)
      throw vtkstd::runtime_error("Missing value.");
    }

  // Ensure we loaded enough values ...
  if(value_count != non_null_size)
    throw vtkstd::runtime_error("Stream doesn't contain enough values.");

  array->Register(0);
  return array;
}

template<typename ValueT>
vtkDenseArray<ValueT>* ReadDenseArrayAscii(istream& stream)
{
  // Create the array ...
  vtkSmartPointer<vtkDenseArray<ValueT> > array = vtkSmartPointer<vtkDenseArray<ValueT> >::New();

  // Read the file header ...
  vtkArrayExtents extents;
  vtkArrayExtents::SizeT non_null_size = 0;
  ReadHeader(stream, extents, non_null_size, array);

  if(non_null_size != extents.GetSize())
    throw vtkstd::runtime_error("Incorrect number of values for a dense array.");

  // Read the file contents ...
  ValueT value;
  vtkArray::SizeT n = 0;
  vtkArrayCoordinates coordinates;
  for(ExtractValue(stream, value); stream; ExtractValue(stream, value), ++n)
    {
    if(n + 1 > non_null_size)
      throw vtkstd::runtime_error("Stream contains too many values.");

    extents.GetRightToLeftCoordinatesN(n, coordinates);
    array->SetValue(coordinates, value);
    }

  if(n != non_null_size)
    throw vtkstd::runtime_error("Stream doesn't contain enough values.");

  array->Register(0);
  return array;
}

} // End anonymous namespace

vtkArrayReader::vtkArrayReader() :
  FileName(0)
{
  this->SetNumberOfInputPorts(0);
}

vtkArrayReader::~vtkArrayReader()
{
  this->SetFileName(0);
}

void vtkArrayReader::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
  os << indent << "FileName: "
     << (this->FileName ? this->FileName : "(none)") << endl;
}

int vtkArrayReader::RequestData(
  vtkInformation*,
  vtkInformationVector**,
  vtkInformationVector* outputVector)
{
  try
    {
    if(!this->FileName)
      throw vtkstd::runtime_error("FileName not set.");

    ifstream file(this->FileName);

    vtkArray* const array = this->Read(file);
    if(!array)
      throw vtkstd::runtime_error("Error reading array.");

    vtkArrayData* const array_data = vtkArrayData::GetData(outputVector);
    array_data->ClearArrays();
    array_data->AddArray(array);
    array->Delete();

    return 1;
    }
  catch(vtkstd::exception& e)
    {
    vtkErrorMacro(<< e.what());
    }

  return 0;
}

vtkArray* vtkArrayReader::Read(istream& stream)
{
  try
    {
    // Read enough of the file header to identify the type ...
    vtkstd::string header_string;
    vtkstd::getline(stream, header_string);
    vtkstd::istringstream header_buffer(header_string);

    vtkstd::string header_magic;
    vtkstd::string header_type;
    header_buffer >> header_magic >> header_type;

    // Read input file type, binary or ascii
    vtkstd::string header_file_string;
    vtkstd::string header_file_type;
    vtkstd::getline(stream, header_file_string);
    vtkstd::istringstream header_file_type_buffer(header_file_string);
    header_file_type_buffer >> header_file_type;

    bool read_binary = false;
    if(header_file_type == "binary")
      {
      read_binary = true;
      }
    else if(header_file_type != "ascii")
      {
      throw vtkstd::runtime_error("Unknown file type: " + header_file_type);
      }

    if(header_magic == "vtk-sparse-array")
      {
      if(header_type == "integer")
        {
        return (read_binary ? ReadSparseArrayBinary<vtkIdType>(stream) : ReadSparseArrayAscii<vtkIdType>(stream));
        }
      else if(header_type == "double")
        {
        return (read_binary ? ReadSparseArrayBinary<double>(stream) : ReadSparseArrayAscii<double>(stream));
        }
      else if(header_type == "string")
        {
        return (read_binary ? ReadSparseArrayBinary<vtkStdString>(stream) : ReadSparseArrayAscii<vtkStdString>(stream));
        }
      else if(header_type == "unicode-string")
        {
        return (read_binary ? ReadSparseArrayBinary<vtkUnicodeString>(stream) : ReadSparseArrayAscii<vtkUnicodeString>(stream));
        }
      else
        {
        throw vtkstd::runtime_error("Unknown array type: " + header_type);
        }
      }
    else if(header_magic == "vtk-dense-array")
      {
      if(header_type == "integer")
        {
        return (read_binary ? ReadDenseArrayBinary<vtkIdType>(stream) : ReadDenseArrayAscii<vtkIdType>(stream));
        }
      else if(header_type == "double")
        {
        return (read_binary ? ReadDenseArrayBinary<double>(stream) : ReadDenseArrayAscii<double>(stream));
        }
      else if(header_type == "string")
        {
        return (read_binary ? ReadDenseArrayBinary<vtkStdString>(stream) : ReadDenseArrayAscii<vtkStdString>(stream));
        }
      else if(header_type == "unicode-string")
        {
        return (read_binary ? ReadDenseArrayBinary<vtkUnicodeString>(stream) : ReadDenseArrayAscii<vtkUnicodeString>(stream));
        }
      else
        {
        throw vtkstd::runtime_error("Unknown array type: " + header_type);
        }
      }
    else
      {
      throw vtkstd::runtime_error("Unknown file type: " + header_magic);
      }
    }
  catch(vtkstd::exception& e)
    {
    vtkGenericWarningMacro(<< e.what());
    }

  return 0;
}