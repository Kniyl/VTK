/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkMergePoints.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkMergePoints.h"

#include "vtkDataArray.h"
#include "vtkIdList.h"
#include "vtkObjectFactory.h"
#include "vtkPoints.h"
#include "vtkFloatArray.h"

vtkCxxRevisionMacro(vtkMergePoints, "1.40");
vtkStandardNewMacro(vtkMergePoints);

// Determine whether point given by x[3] has been inserted into points list.
// Return id of previously inserted point if this is true, otherwise return
// -1.
vtkIdType vtkMergePoints::IsInsertedPoint(const double x[3])
{
  int i, ijk0, ijk1, ijk2;
  vtkIdType idx;
  vtkIdList *bucket;
  //
  //  Locate bucket that point is in.
  //
  ijk0 = (int) ((double) ((x[0] - this->Bounds[0]) / 
           (this->Bounds[1] - this->Bounds[0])) * (this->Divisions[0] - 1));
  ijk1 = (int) ((double) ((x[1] - this->Bounds[2]) / 
           (this->Bounds[3] - this->Bounds[2])) * (this->Divisions[1] - 1));
  ijk2 = (int) ((double) ((x[2] - this->Bounds[4]) / 
           (this->Bounds[5] - this->Bounds[4])) * (this->Divisions[2] - 1));


  idx = ijk0 + ijk1*this->Divisions[0] + 
        ijk2*this->Divisions[0]*this->Divisions[1];

  bucket = this->HashTable[idx];

  if ( ! bucket )
    {
    return -1;
    }
  else // see whether we've got duplicate point
    {
    //
    // Check the list of points in that bucket.
    //
    vtkIdType ptId;
    int nbOfIds = bucket->GetNumberOfIds ();

    // For efficiency reasons, we break the data abstraction for points
    // and ids (we are assuming and vtkIdList
    // is storing ints).
    vtkDataArray *dataArray = this->Points->GetData();
    vtkIdType *idArray = bucket->GetPointer(0);
    if (dataArray->GetDataType() == VTK_FLOAT)
      {
      float f[3];
      f[0] = static_cast<float>(x[0]);
      f[1] = static_cast<float>(x[1]);
      f[2] = static_cast<float>(x[2]);
      vtkFloatArray *floatArray = static_cast<vtkFloatArray *>(dataArray);
      float *pt;
      for (i=0; i < nbOfIds; i++) 
        {
        ptId = idArray[i];
        pt = floatArray->GetPointer(0) + 3*ptId;
        if ( f[0] == pt[0] && f[1] == pt[1] && f[2] == pt[2] )
          {
          return ptId;
          }
        }
      }
    else
      {
      // Using the double interface
      double *pt;
      for (i=0; i < nbOfIds; i++) 
        {
        ptId = idArray[i];
        pt = dataArray->GetTuple(ptId);
        if ( x[0] == pt[0] && x[1] == pt[1] && x[2] == pt[2] )
          {
          return ptId;
          }
        }
      }
    }

  return -1;
}


int vtkMergePoints::InsertUniquePoint(const double x[3], vtkIdType &id)
{
  int i, ijk0, ijk1, ijk2;
  vtkIdType idx;
  vtkIdList *bucket;

  //
  //  Locate bucket that point is in.
  //
  ijk0 = (int) ((double) ((x[0] - this->Bounds[0]) / 
           (this->Bounds[1] - this->Bounds[0])) * (this->Divisions[0] - 1));
  ijk1 = (int) ((double) ((x[1] - this->Bounds[2]) / 
           (this->Bounds[3] - this->Bounds[2])) * (this->Divisions[1] - 1));
  ijk2 = (int) ((double) ((x[2] - this->Bounds[4]) / 
           (this->Bounds[5] - this->Bounds[4])) * (this->Divisions[2] - 1));

  idx = ijk0 + ijk1*this->Divisions[0] + 
        ijk2*this->Divisions[0]*this->Divisions[1];

  bucket = this->HashTable[idx];

  if (bucket) // see whether we've got duplicate point
    {
    //
    // Check the list of points in that bucket.
    //
    vtkIdType ptId;
    int nbOfIds = bucket->GetNumberOfIds ();

    // For efficiency reasons, we break the data abstraction for points
    // and ids (we are assuming vtkPoints stores a vtkIdList
    // is storing ints).
    vtkDataArray *dataArray = this->Points->GetData();
    vtkIdType *idArray = bucket->GetPointer(0);
    
    if (dataArray->GetDataType() == VTK_FLOAT)
      {
      float f[3];
      f[0] = static_cast<float>(x[0]);
      f[1] = static_cast<float>(x[1]);
      f[2] = static_cast<float>(x[2]);
      vtkFloatArray *floatArray = static_cast<vtkFloatArray *>(dataArray);
      float *pt;
      for (i=0; i < nbOfIds; i++) 
        {
        ptId = idArray[i];
        pt = floatArray->GetPointer(0) + 3*ptId;
        if ( f[0] == pt[0] && f[1] == pt[1] && f[2] == pt[2] )
          {
          // point is already in the list, return 0 and set the id parameter
          id = ptId;
          return 0;
          }
        }
      }
    else
      {
      // Using the double interface
      double *pt;
      for (i=0; i < nbOfIds; i++) 
        {
        ptId = idArray[i];
        pt = dataArray->GetTuple(ptId);
        if ( x[0] == pt[0] && x[1] == pt[1] && x[2] == pt[2] )
          {
          // point is already in the list, return 0 and set the id parameter
          id = ptId;
          return 0;
          }
        }
      }
    }
  else
    {
    // create a bucket point list and insert the point
    bucket = vtkIdList::New();
    bucket->Allocate(this->NumberOfPointsPerBucket/2,
                     this->NumberOfPointsPerBucket/3);
    this->HashTable[idx] = bucket;
    }

  // point has to be added
  bucket->InsertNextId(this->InsertionPointId);
  this->Points->InsertPoint(this->InsertionPointId,x);
  id = this->InsertionPointId++;

  return 1;
}

//----------------------------------------------------------------------------
void vtkMergePoints::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}
