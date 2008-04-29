/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkMutableUndirectedGraph.h

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
// .NAME vtkMutableUndirectedGraph - An editable undirected graph.
//
// .SECTION Description
// vtkMutableUndirectedGraph is an undirected graph with additional functions
// for adding vertices and edges. ShallowCopy(), DeepCopy(), CheckedShallowCopy(),
// and CheckedDeepCopy() will succeed when the argument is a vtkUndirectedGraph
// or vtkMutableUndirectedGraph.
//
// .SECTION See Also
// vtkUndirectedGraph vtkGraph

#ifndef __vtkMutableUndirectedGraph_h
#define __vtkMutableUndirectedGraph_h

#include "vtkUndirectedGraph.h"

class vtkEdgeListIterator;
class vtkGraphEdge;

class VTK_FILTERING_EXPORT vtkMutableUndirectedGraph : public vtkUndirectedGraph
{
public:
  static vtkMutableUndirectedGraph *New();
  vtkTypeRevisionMacro(vtkMutableUndirectedGraph, vtkUndirectedGraph);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Adds a new vertex to the graph and returns the id of that vertex.
  vtkIdType AddVertex();

  // Description:
  // Adds a vertex with the given name to the graph (if a vertex of
  // that name does not already exist) and returns the id the vertex
  // with that name.
  vtkIdType AddVertex(const vtkVariant& name);

  // Description:
  // Adds a vertex with the given name to the graph (if a vertex of
  // that name does not already exist). If vertex is non-NULL, the ID
  // of the vertex with the given name will be written into *vertex.
  // In a distributed graph, passing NULL for the second argument
  // can improve performance when adding non-local edges.
  void AddVertex(const vtkVariant& name, vtkIdType *vertex);

  //BTX
  // Description:
  // Adds a new vertex, with properties, to the graph and returns the id of that vertex.
  vtkIdType AddVertex(vtkVariantArray *variantValueArr);

  // Description:
  // Adds an undirected edge between u and v, and returns
  // a vtkEdgeType structure for that edge.
  // The returned vtkEdgeType indicates a Source and Target,
  // but these are in arbitrary order.
  vtkEdgeType AddEdge(vtkIdType u, vtkIdType v);
  
  // Description:
  // Adds an undirected edge, with properties, between u and v, and returns
  // a vtkEdgeType structure for that edge.
  // The returned vtkEdgeType indicates a Source and Target,
  // but these are in arbitrary order.
  vtkEdgeType AddEdge(vtkIdType u, vtkIdType v, vtkVariantArray *variantValueArr);

  // Description:
  // Adds a directed edge from u to v to the graph and returns
  // a vtkEdgeType structure for that edge. u is the name of a
  // vertex, which will be automatically added if it does not
  // already exist.
  vtkEdgeType AddEdge(const vtkVariant& uName, vtkIdType v);

  // Description:
  // Adds a directed edge from u to v to the graph and returns
  // a vtkEdgeType structure for that edge. v is the name of a
  // vertex, which will be automatically added if it does not
  // already exist.
  vtkEdgeType AddEdge(vtkIdType u, const vtkVariant& vName);

  // Description:
  // Adds a directed edge from u to v to the graph and returns
  // a vtkEdgeType structure for that edge. u and v are the names
  // of vertices, which will be automatically added if they do
  // not already exist.
  vtkEdgeType AddEdge(const vtkVariant& uName, const vtkVariant& vName);

  // Description:
  // Adds an undirected edge from u to v to the graph. If non-null, edge
  // will receive the newly-constructed edge. For distributed graphs, 
  // passing NULL for edge can improve performance when adding non-local
  // edges.
  void AddEdge(vtkIdType u, vtkIdType v, vtkEdgeType *edge);

  // Description:
  // Adds a directed edge from u to v to the graph. If non-null, edge
  // will receive the newly-constructed edge. For distributed graphs, 
  // passing NULL for edge can improve performance when adding non-local
  // edges. u is the name of a vertex, which will be automatically
  // added if it does not already exist.
  void AddEdge(const vtkVariant& uName, vtkIdType v, vtkEdgeType *edge);

  // Description:
  // Adds a directed edge from u to v to the graph. If non-null, edge
  // will receive the newly-constructed edge. For distributed graphs, 
  // passing NULL for edge can improve performance when adding non-local
  // edges. v is the name of a vertex, which will be automatically
  // added if it does not already exist.
  void AddEdge(vtkIdType u, const vtkVariant& vName, vtkEdgeType *edge);

  // Description:
  // Adds a directed edge from u to v to the graph. If non-null, edge
  // will receive the newly-constructed edge. For distributed graphs, 
  // passing NULL for edge can improve performance when adding non-local
  // edges. u and v are the names of vertices, which will be automatically
  // added if they do not already exist.
  void AddEdge(const vtkVariant& uName, const vtkVariant& vName, 
               vtkEdgeType *edge);
  
  // Description:
  // Adds an udirected edge, with properties, from u to v to the graph. If non-null, edge
  // will receive the newly-constructed edge. For distributed graphs, 
  // passing NULL for edge can improve performance when adding non-local
  // edges.
  void AddEdge(vtkIdType u, vtkIdType v, vtkEdgeType *edge, vtkVariantArray *variantValueArr);
  //ETX

  // Description:
  // Version of AddEdge that returns a heavyweight vtkGraphEdge
  // for use with wrappers.
  // The graph owns the reference of the edge and will replace
  // its contents on the next call to AddGraphEdge.
  vtkGraphEdge *AddGraphEdge(vtkIdType u, vtkIdType v);

protected:
  vtkMutableUndirectedGraph();
  ~vtkMutableUndirectedGraph();

  // Description:
  // Graph edge that is reused of AddGraphEdge calls.
  vtkGraphEdge *GraphEdge;

private:
  vtkMutableUndirectedGraph(const vtkMutableUndirectedGraph&);  // Not implemented.
  void operator=(const vtkMutableUndirectedGraph&);  // Not implemented.
};

#endif
