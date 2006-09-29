/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkTemporalDataSetCache.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkTemporalDataSetCache - cache time steps
// .SECTION Description
// vtkTemporalDataSetCache cache time step requests of a temporal dataset,
// when cached data is requested it is returned using a shallow copy.


#ifndef __vtkTemporalDataSetCache_h
#define __vtkTemporalDataSetCache_h

#include "vtkTemporalDataSetAlgorithm.h"

#include <vtkstd/map>

class VTK_HYBRID_EXPORT vtkTemporalDataSetCache : public vtkTemporalDataSetAlgorithm
{
public:
  static vtkTemporalDataSetCache *New();
  vtkTypeRevisionMacro(vtkTemporalDataSetCache, vtkTemporalDataSetAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // This is the maximum number of time steps that can be retained in memory.
  // it defaults to 10.
  void SetCacheSize(int size);
  vtkGetMacro(CacheSize,int);

protected:
  vtkTemporalDataSetCache();
  ~vtkTemporalDataSetCache();

  int CacheSize;

//BTX
  typedef vtkstd::map<double,vtkstd::pair<unsigned long,vtkDataObject *> >
  CacheType;
  CacheType Cache;
//ETX
  virtual int RequestUpdateExtent (vtkInformation *,
                                   vtkInformationVector **,
                                   vtkInformationVector *);
  
  virtual int RequestData(vtkInformation *,
                          vtkInformationVector **,
                          vtkInformationVector *);

private:
  vtkTemporalDataSetCache(const vtkTemporalDataSetCache&);  // Not implemented.
  void operator=(const vtkTemporalDataSetCache&);  // Not implemented.
};



#endif



