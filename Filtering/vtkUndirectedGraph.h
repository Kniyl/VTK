/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkUndirectedGraph.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
/*----------------------------------------------------------------------------
 Copyright (c) Sandia Corporation
 See Copyright.txt or http://www.paraview.org/HTML/Copyright.html for details.
----------------------------------------------------------------------------*/
// .NAME vtkUndirectedGraph - An undirected graph.
//
// .SECTION Description
// vtkUndirectedGraph is a collection of vertices along with a collection of
// undirected edges (they connect two vertices in no particular order). 
// ShallowCopy(), DeepCopy(), CheckedShallowCopy(), CheckedDeepCopy()
// accept instances of vtkUndirectedGraph and vtkMutableUndirectedGraph.
// GetOutEdges(v, it) and GetInEdges(v, it) return the same list of edges, 
// which is the list of all edges which have a v as an endpoint.
// GetInDegree(v), GetOutDegree(v) and GetDegree(v) all return the full
// degree of vertex v.
//
// vtkUndirectedGraph is read-only. To create an undirected graph,
// use an instance of vtkMutableUndirectedGraph, then you may set the
// structure to a vtkUndirectedGraph using ShallowCopy().
//
// .SECTION See Also
// vtkGraph vtkMutableUndirectedGraph

#ifndef __vtkUndirectedGraph_h
#define __vtkUndirectedGraph_h

#include "vtkGraph.h"

class VTK_FILTERING_EXPORT vtkUndirectedGraph : public vtkGraph
{
public:
  static vtkUndirectedGraph *New();
  vtkTypeRevisionMacro(vtkUndirectedGraph, vtkGraph);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Return what type of dataset this is.
  virtual int GetDataObjectType() {return VTK_UNDIRECTED_GRAPH;}

  // Description:
  // Returns the full degree of the vertex.
  virtual vtkIdType GetInDegree(vtkIdType v);

  //BTX
  // Description:
  // Retrieve a graph from an information vector.
  static vtkUndirectedGraph *GetData(vtkInformation *info);
  static vtkUndirectedGraph *GetData(vtkInformationVector *v, int i=0);
  //ETX

protected:
  vtkUndirectedGraph();
  ~vtkUndirectedGraph();

  // Description:
  // Check the structure, and accept it if it is a valid
  // undirected graph.
  virtual bool IsStructureValid(vtkGraph *g);

  //BTX
  // Description:
  // For iterators, returns the same edge list as GetOutEdges().
  virtual void GetInEdges(vtkIdType v, const vtkInEdgeType *& edges, vtkIdType & nedges);
  //ETX

private:
  vtkUndirectedGraph(const vtkUndirectedGraph&);  // Not implemented.
  void operator=(const vtkUndirectedGraph&);  // Not implemented.
};

#endif
