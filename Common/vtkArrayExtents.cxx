/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkArrayExtents.cxx
  
-------------------------------------------------------------------------
  Copyright 2008 Sandia Corporation.
  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
  the U.S. Government retains certain rights in this software.
-------------------------------------------------------------------------

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkArrayCoordinates.h"
#include "vtkArrayExtents.h"

#include <vtksys/stl/functional>
#include <vtksys/stl/numeric>

vtkArrayExtents::vtkArrayExtents()
{
}

vtkArrayExtents::vtkArrayExtents(const vtkIdType i) :
  Storage(1)
{
  this->Storage[0] = vtkArrayRange(0, i);
}

vtkArrayExtents::vtkArrayExtents(const vtkArrayRange& i) :
  Storage(1)
{
  this->Storage[0] = i;
}

vtkArrayExtents::vtkArrayExtents(const vtkIdType i, const vtkIdType j) :
  Storage(2)
{
  this->Storage[0] = vtkArrayRange(0, i);
  this->Storage[1] = vtkArrayRange(0, j);
}

vtkArrayExtents::vtkArrayExtents(const vtkArrayRange& i, const vtkArrayRange& j) :
  Storage(2)
{
  this->Storage[0] = i;
  this->Storage[1] = j;
}

vtkArrayExtents::vtkArrayExtents(const vtkIdType i, const vtkIdType j, const vtkIdType k) :
  Storage(3)
{
  this->Storage[0] = vtkArrayRange(0, i);
  this->Storage[1] = vtkArrayRange(0, j);
  this->Storage[2] = vtkArrayRange(0, k);
}

vtkArrayExtents::vtkArrayExtents(const vtkArrayRange& i, const vtkArrayRange& j, const vtkArrayRange& k) :
  Storage(3)
{
  this->Storage[0] = i;
  this->Storage[1] = j;
  this->Storage[2] = k;
}

const vtkArrayExtents vtkArrayExtents::Uniform(vtkIdType n, vtkIdType m)
{
  vtkArrayExtents result;
  // IA64 HP-UX doesnt seem to have the vector<T> vector1(n, value) 
  // overload nor the assign(n, value) method, so we use the single 
  // arguement constructor and initialize the values manually.
  // result.Storage = vtksys_stl::vector<vtkIdType>(n, m);

  result.Storage = vtksys_stl::vector<vtkArrayRange>(n);
  for(vtkIdType i = 0; i < n; i++)
    {
    result.Storage[i] = vtkArrayRange(0, m);
    }
  return result;
}

void vtkArrayExtents::Append(const vtkArrayRange& extent)
{
  this->Storage.push_back(extent);
}

vtkIdType vtkArrayExtents::GetDimensions() const
{
  return this->Storage.size();
}

vtkIdType vtkArrayExtents::GetSize() const
{
  if(this->Storage.empty())
    return 0;

  vtkIdType size = 1;
  for(vtkIdType i = 0; i != this->Storage.size(); ++i)
    size *= this->Storage[i].GetSize();

  return size;
}

void vtkArrayExtents::SetDimensions(vtkIdType dimensions)
{
  this->Storage.assign(dimensions, vtkArrayRange());
}

vtkArrayRange& vtkArrayExtents::operator[](vtkIdType i)
{
  return this->Storage[i];
}

const vtkArrayRange& vtkArrayExtents::operator[](vtkIdType i) const
{
  return this->Storage[i];
}

bool vtkArrayExtents::operator==(const vtkArrayExtents& rhs) const
{
  return this->Storage == rhs.Storage;
}

bool vtkArrayExtents::operator!=(const vtkArrayExtents& rhs) const
{
  return !(*this == rhs);
}

bool vtkArrayExtents::ZeroBased() const
{
  for(vtkIdType i = 0; i != this->GetDimensions(); ++i)
    {
    if(this->Storage[i].GetBegin() != 0)
      return false;
    }

  return true;
}

bool vtkArrayExtents::SameShape(const vtkArrayExtents& rhs) const
{
  if(this->GetDimensions() != rhs.GetDimensions())
    return false;

  for(vtkIdType i = 0; i != this->GetDimensions(); ++i)
    {
    if(this->Storage[i].GetSize() != rhs.Storage[i].GetSize())
      return false;
    }

  return true;
}

void vtkArrayExtents::GetLeftToRightCoordinatesN(vtkIdType n, vtkArrayCoordinates& coordinates) const
{
  coordinates.SetDimensions(this->GetDimensions());

  vtkIdType divisor = 1;
  for(vtkIdType i = 0; i < this->GetDimensions(); ++i)
    {
    coordinates[i] = ((n / divisor) % this->Storage[i].GetSize()) + this->Storage[i].GetBegin();
    divisor *= this->Storage[i].GetSize();
    }
}

void vtkArrayExtents::GetRightToLeftCoordinatesN(vtkIdType n, vtkArrayCoordinates& coordinates) const
{
  coordinates.SetDimensions(this->GetDimensions());

  vtkIdType divisor = 1;
  for(vtkIdType i = this->GetDimensions() - 1; i >= 0; --i)
    {
    coordinates[i] = ((n / divisor) % this->Storage[i].GetSize()) + this->Storage[i].GetBegin();
    divisor *= this->Storage[i].GetSize();
    }
}

bool vtkArrayExtents::Contains(const vtkArrayExtents& other) const
{
  if(this->GetDimensions() != other.GetDimensions())
    return false;

  for(vtkIdType i = 0; i != this->GetDimensions(); ++i)
    {
    if(!this->Storage[i].Contains(other[i]))
      return false;
    }

  return true;
}

bool vtkArrayExtents::Contains(const vtkArrayCoordinates& coordinates) const
{
  if(coordinates.GetDimensions() != this->GetDimensions())
    return false;

  for(vtkIdType i = 0; i != this->GetDimensions(); ++i)
    {
    if(!this->Storage[i].Contains(coordinates[i]))
      return false;
    }

  return true;
}

ostream& operator<<(ostream& stream, const vtkArrayExtents& rhs)
{
  for(unsigned int i = 0; i != rhs.Storage.size(); ++i)
    {
    if(i)
      stream << "x";
    stream << "[" << rhs.Storage[i].GetBegin() << "," << rhs.Storage[i].GetEnd() << ")";
    }

  return stream;
}

