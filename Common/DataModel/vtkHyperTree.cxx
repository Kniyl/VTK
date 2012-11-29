/*=========================================================================

Program:   Visualization Toolkit
Module:    vtkHyperTreeGrid.cxx

Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
All rights reserved.
See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkHyperTree.h"

#include "vtkHyperTreeCursor.h"
#include "vtkObjectFactory.h"

#include <deque>
#include <vector>

#include <assert.h>

// Description:
// The template value N describes the number of children to binary and
// ternary trees.
template<int N> class vtkCompactHyperTree; // : public vtkHyperTree
template<int N> class vtkCompactHyperTreeNode;
template<int N> class vtkCompactHyperTreeCursor : public vtkHyperTreeCursor
{
public:
  //---------------------------------------------------------------------------
  vtkTypeMacro(vtkCompactHyperTreeCursor<N>,vtkHyperTreeCursor);

  static vtkCompactHyperTreeCursor<N>* New()
  {
    vtkObject* o = vtkObjectFactory::CreateInstance( "vtkCompactHyperTreeCursor<N>" );

    if( o )
      {
      return static_cast<vtkCompactHyperTreeCursor<N> *>( o );
      }
    else
      {
      return new vtkCompactHyperTreeCursor<N>;
      }
  }

  //---------------------------------------------------------------------------
  // Initialization
  virtual void SetTree( vtkCompactHyperTree<N>* tree )
  {
    this->Tree = tree;
  }

  //---------------------------------------------------------------------------
  int GetLeafId()
  {
    assert( "pre: is_leaf" && IsLeaf() );
    return this->Index;
  }

  //---------------------------------------------------------------------------
  virtual bool IsLeaf()
  {
    return this->Leaf;
  }

  //---------------------------------------------------------------------------
  virtual bool IsTerminalNode()
  {
    bool result = ! this->Leaf;
    if( result )
      {
      vtkCompactHyperTreeNode<N>* node = this->Tree->GetNode( this->Index );
      result = node->IsTerminalNode();
      }
    // A=>B: notA or B
    assert( "post: compatible" && ( ! result || ! this->Leaf) );
    return result;
  }

  //---------------------------------------------------------------------------
  virtual bool IsRoot()
  {
    return ( ! this->Leaf && this->Index == 1 )
      || ( this->Leaf && ! this->Index && this->Tree->GetLeafParentSize() == 1 );
  }

  //---------------------------------------------------------------------------
  virtual int GetCurrentLevel()
  {
    int result = this->GetChildHistorySize();
    assert( "post: positive_result" && result >= 0 );
    return result;
  }

  //---------------------------------------------------------------------------
  // Description:
  // Return the child number of the current node relative to its parent.
  // \pre not_root: !IsRoot().
  // \post valid_range: result >= 0 && result<GetNumberOfChildren()
  virtual int GetChildIndex()
  {
    assert( "post: valid_range" && this->ChildIndex >= 0 && this->ChildIndex<GetNumberOfChildren() );
    return this->ChildIndex;
  }

  //---------------------------------------------------------------------------
  // Cursor movement.
  // \pre can be root
  // \post is_root: IsRoot()
  virtual void ToRoot()
  {
    this->ChildHistory.clear();
    this->Leaf = ( this->Tree->GetLeafParentSize() == 1 );
    if ( this->Leaf )
      {
      this->Index = 0;
      }
    else
      {
      this->Index = 1;
      }
    this->ChildIndex = 0;
    
    for ( unsigned int i = 0; i < this->Dimension; ++ i )
      {
      this->Indices[i] = 0;
      }
  }

  //---------------------------------------------------------------------------
  // \pre not_root: !IsRoot()
  virtual void ToParent()
  {
    assert( "pre: not_root" && !IsRoot() );
    if( this->Leaf)
      {
      this->Index = this->Tree->GetLeafParent( this->Index);
      }
    else
      {
      this->Index = this->Tree->GetNode( this->Index)->GetParent();
      }
    this->Leaf = 0;
    this->ChildIndex=this->ChildHistory.back(); // top()
    this->ChildHistory.pop_back();

    for ( unsigned int i = 0; i < this->Dimension;  ++ i )
      {
      this->Indices[i]=( this->Indices[i] ) / this->Tree->GetBranchFactor();
      }
  }

  //---------------------------------------------------------------------------
  // \pre not_leaf: !IsLeaf()
  // \pre valid_child: child >= 0 && child<this->GetNumberOfChildren()
  virtual void ToChild(int child)
  {
    assert( "pre: not_leaf" && !IsLeaf() );
    assert( "pre: valid_child" && child >= 0 && child<this->GetNumberOfChildren() );

    vtkCompactHyperTreeNode<N> *node=this->Tree->GetNode( this->Index);
    this->ChildHistory.push_back( this->ChildIndex );
    this->ChildIndex=child;
    this->Index=node->GetChild( child );
    this->Leaf=node->IsChildLeaf( child );

    int tmpChild = child;
    int tmp;
    int branchFactor = this->Tree->GetBranchFactor();
    for ( unsigned int i = 0; i < this->Dimension; ++ i )
      { 
      // Effectively convert child to base 2/3 (branch factor)
      tmp = tmpChild;
      tmpChild /= branchFactor;
      int index=tmp-(branchFactor*tmpChild); // Remainder (mod)
      assert( "check: mod 3 value" && index >= 0 && index<branchFactor);
      this->Indices[i]=(( this->Indices[i])*branchFactor)+index;
      }
  }

  //---------------------------------------------------------------------------
  // Description:
  // Move the cursor to the same node pointed by `other'.
  // \pre other_exists: other != 0
  // \pre same_hyperTree: this->SameTree( other )
  // \post equal: this->IsEqual( other )
  virtual void ToSameNode( vtkHyperTreeCursor* other )
  {
    assert( "pre: other_exists" && other != 0 );
    assert( "pre: same_hyperTree" && this->SameTree( other ) );

    vtkCompactHyperTreeCursor<N> *o = static_cast<vtkCompactHyperTreeCursor<N> *>( other );

    this->Index = o->Index;
    this->ChildIndex = o->ChildIndex;
    this->Leaf = o->Leaf;
    this->ChildHistory = o->ChildHistory; // use assignment operator
      
    for( unsigned int i = 0; i < this->Dimension; ++ i )
      {
      this->Indices[i] = o->Indices[i];
      }
    assert( "post: equal" && this->IsEqual(other) );
  }

  //--------------------------------------------------------------------------
  // Description:
  // Is `this' equal to `other'?
  // \pre other_exists: other != 0
  // \pre same_hyperTree: this->SameTree(other);
  virtual int IsEqual( vtkHyperTreeCursor* other )
  {
    assert( "pre: other_exists" && other != 0 );
    assert( "pre: same_hyperTree" && this->SameTree(other) );

    vtkCompactHyperTreeCursor<N> *o=static_cast<vtkCompactHyperTreeCursor<N> *>(other);

    int result = this->Index==o->Index && this->ChildIndex==o->ChildIndex
      && this->Leaf==o->Leaf && this->ChildHistory==o->ChildHistory;

    for( unsigned int i = 0; result && i < this->Dimension; ++ i )
      {
      result = this->Indices[i] == o->Indices[i];
      ++ i;
      }
    return result;
  }

  //--------------------------------------------------------------------------
  // Description:
  // Create a copy of `this'.
  // \post results_exists:result != 0
  // \post same_tree: result->SameTree( this )
  virtual vtkHyperTreeCursor* Clone()
  {
    vtkCompactHyperTreeCursor<N>* result = this->NewInstance();
    result->Tree = this->Tree;
    assert( "post: results_exists" && result != 0 );
    assert( "post: same_tree" && result->SameTree( this ) );
    return result;
  }

  //---------------------------------------------------------------------------
  // Description:
  // Are `this' and `other' pointing on the same hyperTree?
  // \pre other_exists: other != 0
  virtual int SameTree( vtkHyperTreeCursor* other )
  {
    assert( "pre: other_exists" && other != 0 );
    vtkCompactHyperTreeCursor<N> *o=vtkCompactHyperTreeCursor<N>::SafeDownCast( other );
    int result = o != 0;
    if(result)
      {
      result = this->Tree==o->Tree;
      }
    return result;
  }

  //---------------------------------------------------------------------------
  // Description:
  // Return the index in dimension `d', as if the node was a cell of a
  // uniform grid of 1<<GetCurrentLevel() cells in each dimension.
  // \pre valid_range: d >= 0 && d<GetDimension()
  // \post valid_result: result >= 0 && result<(1<<GetCurrentLevel() )
  virtual int GetIndex(int d)
  {
    assert( "pre: valid_range" &&  d >= 0 && d<this->Dimension );
    int result = this->Indices[d];
    return result;
  }

  //---------------------------------------------------------------------------
  // Description:
  // Return the number of children for each node of the tree.
  // \post positive_number: result>0
  virtual int GetNumberOfChildren()
  {
    return N;
  }

  //---------------------------------------------------------------------------
  // Description:
  // Return the dimension of the tree.
  // \post positive_result: result >= 0
  virtual int GetDimension()
  {
    assert( "post: positive_result " && this->Dimension>0 );
    assert( "post: up_to_3 " && this->Dimension<=3 ); // and then
    return this->Dimension;
  }

  //---------------------------------------------------------------------------
  // Description:
  // Move to the node described by its indices in each dimension and
  // at a given level. If there is actually a node or a leaf at this
  // location, Found() returns true. Otherwise, Found() returns false and the
  // cursor moves to the closest parent of the query. It can be the root in the
  // worst case.
  // \pre indices_exists: indices != 0
  // \pre valid_size: sizeof(indices)==GetDimension()
  // \pre valid_level: level >= 0
  virtual void MoveToNode(int* indices,
                          int level)
  {
    assert( "pre: indices_exists" && indices != 0 );
    assert( "pre: valid_level" && level >= 0 );

    this->ToRoot();
    int currentLevel = 0;

    int child;
    int tmpIndices[3];

    // Convert to base 2 / 3 starting with most significant digit.
    int mask;
    tmpIndices[0] = indices[0];
    tmpIndices[1] = indices[1];
    tmpIndices[2] = indices[2];
    int i = 0;
    mask = 1;
    while ( ++ i < level )
      {
      mask *= this->Tree->GetBranchFactor();
      }

    while( !this->IsLeaf() && currentLevel < level )
      {
      // Compute the child index
      i = this->Dimension - 1;
      child = 0;
      while ( i >= 0 )
        {
        int digit = tmpIndices[i] / mask;
        tmpIndices[i] -= digit*mask;
        child *= child * this->Tree->GetBranchFactor() + digit;
        -- i;
        }
      this->ToChild( child );
      ++ currentLevel;
      mask /= this->Tree->GetBranchFactor();
      }
    this->IsFound = ( currentLevel == level );
  }

  //---------------------------------------------------------------------------
  // Description
  // Did the last call to MoveToNode succeed?
  virtual int Found()
  {
    return this->IsFound;
  }

  //---------------------------------------------------------------------------
  // Description:
  // Public only for vtkCompactHyperTree.
  void SetIsLeaf( bool value )
  {
    this->Leaf = value;
  }

  //---------------------------------------------------------------------------
  // Description:
  // Public only for vtkCompactHyperTree.
  void SetChildIndex(int childIndex )
  {
    assert( "pre: valid_range" && childIndex >= 0 && childIndex<GetNumberOfChildren() );
    this->ChildIndex = childIndex;
    assert( "post: is_set" && childIndex==GetChildIndex() );
  }

  //---------------------------------------------------------------------------
  // Description:
  // Public only for vtkCompactHyperTree.
  void SetIndex( int index )
  {
    assert( "pre: positive_index" && index >= 0 );
    this->Index = index;
  }

  //---------------------------------------------------------------------------
  // Description:
  // Public only for vtkCompactHyperTree.
  vtkIdType GetChildHistorySize()
  {
    return this->ChildHistory.size();
  }

protected:
  //---------------------------------------------------------------------------
  vtkCompactHyperTreeCursor()
  {
    this->Dimension = 0;
    switch (N)
      {
      case 2:
        this->Dimension = 1;
        break;
      case 3:
        this->Dimension = 1;
        break;
      case 4:
        this->Dimension = 2;
        break;
      case 9:
        this->Dimension = 2;
        break;
      case 8:
        this->Dimension = 3;
        break;
      case 27:
        this->Dimension = 3;
        break;
      default:
        assert( "Bad number of children" && this->Dimension == 0 );
      }
    this->Tree = 0;
    this->Index = 0;
    this->Leaf = 0;
    this->ChildIndex = 0;

    for ( unsigned int i = 0; i < this->Dimension; ++ i )
      {
      this->Indices[i] = 0;
      }
  }

  vtkCompactHyperTree<N> *Tree;
  unsigned char Dimension;

  // Index either in the Nodes or Parents (if leaf)
  int Index;

  // Number of current node as a child
  int ChildIndex;

  int IsFound;
  bool Leaf;

  // A stack, but stack does not have clear()
  std::deque<int> ChildHistory;

  // Index in each dimension of the current node, as if the tree at the current 
  // level were a uniform grid. Default to 3 dimensions, use only those needed
  int Indices[3];

private:
  vtkCompactHyperTreeCursor(const vtkCompactHyperTreeCursor<N> &);  // Not implemented.
  void operator=(const vtkCompactHyperTreeCursor<N> &);    // Not implemented.
};

// We could use a 4 byte int, but the internals are completely hidden.
class vtkHyperTreeLeafFlags
{
public:
  vtkHyperTreeLeafFlags()
  { // Unused bits are set to 1.
    this->Flags[0] = this->Flags[1] = this->Flags[2] = this->Flags[3] = 255;
  }
  // True if all chilren are leaves.
  bool IsTerminal()
  {
    // Unused bits are set to 1.
    return ( this->Flags[0] == 255) && ( this->Flags[1] == 255) && ( this->Flags[2] == 255);
  }
  void SetLeafFlag(int idx, bool val)
  {
    assert( "Valid child idx" && idx >= 0 && idx < 32);
    int i = 0;
    while (idx >= 8)
      {
      ++ i;
      idx-=8;
      }
    unsigned char mask = 1<<idx;
    if ( val )
      {
      this->Flags[i] = this->Flags[i] | mask;
      }
    else
      {
      this->Flags[i] = this->Flags[i] & (mask^255);
      }
  }
  bool GetLeafFlag(int idx )
  {
    assert( "Valid child idx" && idx >= 0 && idx < 32);
    int i = 0;
    while ( idx >= 8 )
      {
      ++ i;
      idx-=8;
      }
    unsigned char mask = 1<<idx;
    return (mask & this->Flags[i]) == mask;
  }
  void PrintSelf(ostream& os, int numChildren)
  {
    assert( "Number of children" && numChildren >= 0 && numChildren < 32);
    int childIdx = 0;
    int byteIdx = 0;
    unsigned char mask = 1;
    while ( childIdx < numChildren )
      {
      os << ((( this->Flags[byteIdx])&mask)==mask);
      ++childIdx;
      if (mask == 128)
        {
        mask = 1;
        ++byteIdx;
        }
      else
        {
        mask<<=1;
        }
      }
    os<<endl;
  }

private:
  unsigned char Flags[4];
};

// Description:
// A node of the Tree which is not a leaf.
// Expected template values: 4, 8, 9, 27.
template<int N> class vtkCompactHyperTreeNode
{
public:
  //---------------------------------------------------------------------------
  // Description:
  // See GetParent().
  void SetParent(int parent)
  {
    assert( "pre: positive_parent" && parent >= 0 );
    this->Parent=parent;
    assert( "post: is_set" && parent==this->GetParent() );
  }

  //---------------------------------------------------------------------------
  // Description:
  // Return the index of the parent node of the current node in the
  // nodes array of the hyperTree.
  int GetParent()
  {
    assert( "post: positive_result" && this->Parent >= 0 );
    return this->Parent;
  }

  //---------------------------------------------------------------------------
  // Description:
  // See GetLeafFlags()
  void SetLeafFlag(int childIdx, bool flag)
  {
    this->LeafFlags.SetLeafFlag(childIdx, flag);
  }

  //---------------------------------------------------------------------------
  // Description
  // Are children all leaves?
  bool IsTerminalNode()
  {
    return this->LeafFlags.IsTerminal();
  }

  //---------------------------------------------------------------------------
  // Description:
  // Is the `i'-th child of the node a leaf ?
  bool IsChildLeaf( int i )
  {
    assert( "pre: valid_range" && i >= 0 && i < N);
    return this->LeafFlags.GetLeafFlag( i );
  }

  //---------------------------------------------------------------------------
  // Description:
  // See GetChild().
  void SetChild( int i, int child )
  {
    assert( "pre: valid_range" && i >= 0 && i < N);
    assert( "pre: positive_child" && child >= 0 );
    this->Children[i] = child;
    assert( "post: is_set" && child==this->GetChild( i ) );
  }

  //---------------------------------------------------------------------------
  // Description:
  // Return the index of of the 'i'-th child. If the result of
  // IsChildLeaf( i ) is true, the index points to an element in the LeafParent
  // and Attribute arrays of the hyperTree class. If not, the index points to
  // an element in the Nodes array of the hyperTree class.
  int GetChild( int i )
  {
    assert( "pre: valid_range" && i >= 0 && i < N);
    assert( "post: positive_result" && this->Children[i] >= 0 );
    return this->Children[i];
  }

  //---------------------------------------------------------------------------
  void PrintSelf(ostream& os, vtkIndent indent)
  {
    os << indent << "Parent=" << this->Parent << endl;
    os << indent << "LeafFlags= ";
    this->LeafFlags.PrintSelf( os, N );

    for( int i = 0; i < N; ++ i )
      {
      os << indent << this->Children[i] << endl;
      }
  }

protected:
  //---------------------------------------------------------------------------
  int Parent; // index
  vtkHyperTreeLeafFlags LeafFlags;
  int Children[N];
};

template<int N> class vtkCompactHyperTree : public vtkHyperTree
{
public:
  vtkTypeMacro(vtkCompactHyperTree<N>,vtkHyperTree);

  //---------------------------------------------------------------------------
  static vtkCompactHyperTree<N>* New()
  {
    vtkObject* o = vtkObjectFactory::CreateInstance( "vtkCompactHyperTree<N>" );
    
    if( o )
      {
      return static_cast<vtkCompactHyperTree<N> *>( o );
      }
    else
      {
      return new vtkCompactHyperTree<N>;
      }
  }

  //---------------------------------------------------------------------------
  // Description:
  // Restore the initial state: only one node and one leaf: the root.
  virtual void Initialize()
  {
    this->Nodes.resize( 1 );
    this->Nodes[0].SetParent( 0 );
    for ( int i = 0; i < N; ++ i )
      {
      // It is assumed that the root is a special node with only one child.
      // The other children flags are irrelevant, but set them as nodes for no good reason.
      this->Nodes[0].SetLeafFlag( i, i == 0 ); // First child is a leaf
      this->Nodes[0].SetChild( i, 0 );
      }
    this->LeafParent.resize( 1 );
    this->LeafParent[0] = 0;
    this->NumberOfLevels = 1;
    this->NumberOfLeavesPerLevel.resize( 1 );
    this->NumberOfLeavesPerLevel[0] = 1;
  }

  //---------------------------------------------------------------------------
  virtual vtkHyperTreeCursor* NewCursor()
  {
    vtkCompactHyperTreeCursor<N>* result = vtkCompactHyperTreeCursor<N>::New();
    result->SetTree( this );
    return result;
  }
  
  //---------------------------------------------------------------------------
  virtual ~vtkCompactHyperTree()
  {
  }
  
  //---------------------------------------------------------------------------
  virtual vtkIdType GetNumberOfLeaves()
  {
    return this->LeafParent.size();
  }
  
  //---------------------------------------------------------------------------
  // Description:
  // Return the erenumber of levels.
  // \post result_greater_or_equal_to_one: result>=1
  virtual vtkIdType GetNumberOfLevels()
  {
    assert( "post: result_greater_or_equal_to_one" && this->NumberOfLevels >= 1);
    return this->NumberOfLevels;
  }
  
  //---------------------------------------------------------------------------
  // Description:
  // Public only for the vtkCompactHyperTreeCursor.
  vtkCompactHyperTreeNode<N>* GetNode( int nodeIdx )
  {
    assert( "pre: valid_range" && nodeIdx >= 0 && nodeIdx<GetNumberOfNodes() );
    return &this->Nodes[nodeIdx];
  }
  
  //---------------------------------------------------------------------------
  // Description:
  // Public only for the vtkCompactHyperTreeCursor.
  // NB: Cursor (index ) appears to be different between nodes and leaves.
  // Different arrays => overlapping indexes.
  // I am changing the name for clarity.
  // This really returns the nodeIdx of the leafs parent.
  int GetLeafParent( int leafIdx )
  {
    assert( "pre: valid_range" && leafIdx >= 0 && leafIdx<this->GetNumberOfLeaves() );
    assert( "post: valid_result" && this->LeafParent[leafIdx] >= 0 && this->LeafParent[leafIdx]<this->GetNumberOfNodes() );
    return this->LeafParent[leafIdx];
  }
  
  //---------------------------------------------------------------------------
  // Description:
  // Public only for the vtkCompactHyperTreeCursor.
  virtual int GetNumberOfNodes()
  {
    assert( "post: not_empty" && this->Nodes.size()>0 );
    return static_cast<int>( this->Nodes.size() );
  }
  
  //---------------------------------------------------------------------------
  // Description:
  // Subdivide node pointed by cursor, only if its a leaf.
  // At the end, cursor points on the node that used to be leaf.
  // \pre leaf_exists: leaf != 0
  // \pre is_a_leaf: leaf->IsLeaf()
  void SubdivideLeaf( vtkHyperTreeCursor* leafCursor )
  {
    assert( "pre: leaf_exists" && leafCursor != 0 );
    assert( "pre: is_a_leaf" && leafCursor->IsLeaf() );

    // We are using a vtkCompactHyperTreeCursor.
    // We know that GetLeafId() return Cursor.
    int leafIndex = leafCursor->GetLeafId();
    vtkCompactHyperTreeCursor<N>* cursor = static_cast<vtkCompactHyperTreeCursor<N> *>(leafCursor);

    // The leaf becomes a node and is not anymore a leaf
    cursor->SetIsLeaf( 0 ); // let the cursor know about that change.
    size_t nodeIndex = this->Nodes.size();
    cursor->SetIndex( static_cast<int>( nodeIndex ) );

    // Nodes get constructed with leaf flags set to 1.
    this->Nodes.resize( nodeIndex + 1 );
    int parentNodeIdx = this->LeafParent[leafIndex];
    this->Nodes[nodeIndex].SetParent( parentNodeIdx );

    // Change the parent: it has one less child as a leaf
    vtkCompactHyperTreeNode<N>* parent = &( this->Nodes[parentNodeIdx] );

    // New nodes index in parents children array.
    int idx = cursor->GetChildIndex();
    assert( "check matching_child" && parent->GetChild( idx ) == leafIndex );
    parent->SetLeafFlag( idx, false );
    parent->SetChild( idx, static_cast<int>( nodeIndex ) );

    // The first new child
    // Recycle the leaf index we are deleting because it became a node.
    // This avoids messy leaf parent array issues.
    this->Nodes[nodeIndex].SetChild( 0, leafIndex );
    this->LeafParent[leafIndex] = static_cast<int>( nodeIndex );

    // The other (N-1) new children.
    size_t nextLeaf = this->LeafParent.size();
    this->LeafParent.resize( nextLeaf + ( N - 1 ) );
    for( int i = 1; i < N; ++ i, ++ nextLeaf )
      {
      this->Nodes[nodeIndex].SetChild( i, static_cast<int>( nextLeaf ) );
      this->LeafParent[nextLeaf] = static_cast<int>( nodeIndex );
      }

    // Update the number of leaves per level.
    int level = cursor->GetChildHistorySize();

    // Remove the subdivided leaf from the number of leaves at its level.
    -- this->NumberOfLeavesPerLevel[level];

    // Add the new leaves to the number of leaves at the next level.
    if( level + 1 == this->NumberOfLevels ) // >=
      {
      // We have a new level.
      ++ this->NumberOfLevels;
      this->NumberOfLeavesPerLevel.resize( this->NumberOfLevels );
      }
    this->NumberOfLeavesPerLevel[level + 1] += N;
  }

  //---------------------------------------------------------------------------
  // NB: Bad interface: This is really GetNumberOfLeaves.
  int GetLeafParentSize()
  {
    return static_cast<int>( this->LeafParent.size() );
  }

  //---------------------------------------------------------------------------
  void PrintSelf( ostream& os, vtkIndent indent )
  {
    this->Superclass::PrintSelf(os,indent);

    os << indent << "Nodes=" << this->Nodes.size() << endl;
    os << indent << "LeafParent=" << this->LeafParent.size() << endl;

    os << indent << "Nodes=" << this->Nodes.size() << endl;
    for ( unsigned int i = 0; i < this->Nodes.size(); ++ i )
      {
      this->Nodes[i].PrintSelf( os, indent );
      }
    os << endl;

    os << indent << "LeafParent="<<this->LeafParent.size() << endl;
    for ( unsigned int i = 0; i < this->LeafParent.size(); ++ i )
      {
      os << this->LeafParent[i] << " ";
      ++ i;
      }
    os << endl;
  }

  //---------------------------------------------------------------------------
  // Description:
  // Return memory used in kilobytes.
  // Ignore the attribute array because its size is added by the data set.
  unsigned int GetActualMemorySize()
  {
    size_t size;
    size = sizeof(int) * this->GetNumberOfLeaves();
    size += sizeof(vtkCompactHyperTreeNode<N>) * this->Nodes.size();
    return static_cast<unsigned int>(size / 1024);
  }

  int GetBranchFactor()
  {
    return this->BranchFactor;
  }

  int GetDimension()
  {
    return this->Dimension;
  }

protected:
  //---------------------------------------------------------------------------
  // Description:
  // Default constructor.
  // The tree as only one node and one leaf: the root.
  vtkCompactHyperTree()
  {
    switch ( N )
      {
      case 2:
        this->BranchFactor = 2;
        this->Dimension = 1;
        break;
      case 3:
        this->BranchFactor = 3;
        this->Dimension = 1;
        break;
      case 4:
        this->BranchFactor = 2;
        this->Dimension = 2;
        break;
      case 8:
        this->BranchFactor = 2;
        this->Dimension = 3;
        break;
      case 9:
        this->BranchFactor = 3;
        this->Dimension = 2;
        break;
      case 27:
        this->BranchFactor = 3;
        this->Dimension = 3;
        break;
      } // switch ( N )

    // Set root node
    this->Nodes.resize( 1 );
    this->Nodes[0].SetParent( 0 );

    // Nodes default to have all children leaf flags equal true
    for ( int i = 0 ; i < N ; ++ i )
      {
      this->Nodes[0].SetLeafFlag( i, i == 0 ); // First child is a leaf
      this->Nodes[0].SetChild( i, 0 );
      }
    this->LeafParent.resize( 1 );
    this->LeafParent[0] = 0;
    this->NumberOfLevels = 1;
    this->NumberOfLeavesPerLevel.resize( 1 );
    this->NumberOfLeavesPerLevel[0] = 1;
  }

  int BranchFactor;
  int Dimension;
  vtkIdType NumberOfLevels;
  std::vector<vtkCompactHyperTreeNode<N> > Nodes;

  // Storage for number of leaves in each level
  std::vector<int> NumberOfLeavesPerLevel;

  // Storage to record the parent of each leaf
  std::vector<int> LeafParent;

private:
  vtkCompactHyperTree(const vtkCompactHyperTree<N> &);  // Not implemented.
  void operator=(const vtkCompactHyperTree<N> &);    // Not implemented.
};

//-----------------------------------------------------------------------------
vtkHyperTree* vtkHyperTree::CreateInstance( int factor, int dimension )
{
  switch ( factor )
    {
    case 2:
      switch( dimension )
        {
        case 3:
          return vtkCompactHyperTree<8>::New();
        case 2:
          return vtkCompactHyperTree<4>::New();
        case 1:
          return vtkCompactHyperTree<2>::New();
        default:
          vtkGenericWarningMacro( "Bad dimension " << dimension );
          return NULL;
        }
    case 3:
      switch( dimension )
        {
        case 3:
          return vtkCompactHyperTree<27>::New();
        case 2:
          return vtkCompactHyperTree<9>::New();
        case 1:
          return vtkCompactHyperTree<3>::New();
        default:
          vtkGenericWarningMacro( "Bad dimension " << dimension );
          return NULL;
        }
    default:
      vtkGenericWarningMacro( "Bad branching factor " << factor);
    }

  return NULL;
}

//-----------------------------------------------------------------------------
void vtkHyperTree::FindChildParameters( int child, int &index, bool& isLeaf )
{
  
  switch ( this->GetDimension() )
    {
    case 3:
      switch ( this->GetBranchFactor() )
        {
        case 2:
          {
          vtkCompactHyperTree<8>* tree;
          tree = static_cast<vtkCompactHyperTree<8>*>( this );
          vtkCompactHyperTreeNode<8>* node = tree->GetNode( index );
          index = node->GetChild( child );
          isLeaf = node->IsChildLeaf( child );
          break;
          }
        case 3:
          {
          vtkCompactHyperTree<27>* tree;
          tree = static_cast<vtkCompactHyperTree<27>*>( this );
          vtkCompactHyperTreeNode<27>* node = tree->GetNode( index );
          index = node->GetChild( child );
          isLeaf = node->IsChildLeaf( child );
          break;
          }
        default:
          vtkGenericWarningMacro( "Bad branching factor " << this->GetBranchFactor() );
          return;
        } // case 3
      break;
    case 2:
      switch ( this->GetBranchFactor() )
        {
        case 2:
          {
          vtkCompactHyperTree<4>* tree;
          tree = static_cast<vtkCompactHyperTree<4>*>( this );
          vtkCompactHyperTreeNode<4>* node = tree->GetNode( index );
          index = node->GetChild( child );
          isLeaf = node->IsChildLeaf( child );
          break;
          }
        case 3:
          {
          vtkCompactHyperTree<9>* tree;
          tree = static_cast<vtkCompactHyperTree<9>*>( this );
          vtkCompactHyperTreeNode<9>* node = tree->GetNode( index );
          index = node->GetChild( child );
          isLeaf = node->IsChildLeaf( child );
          break;
          }
        default:
          vtkGenericWarningMacro( "Bad branching factor " << this->GetBranchFactor() );
          return;
        } // case 2
      break;
    case 1:
      switch ( this->GetBranchFactor() )
        {
        case 2:
          {
          vtkCompactHyperTree<2>* tree;
          tree = static_cast<vtkCompactHyperTree<2>*>( this );
          vtkCompactHyperTreeNode<2>* node = tree->GetNode( index );
          index = node->GetChild( child );
          isLeaf = node->IsChildLeaf( child );
          break;
          }
        case 3:
          {
          vtkCompactHyperTree<3>* tree;
          tree = static_cast<vtkCompactHyperTree<3>*>( this );
          vtkCompactHyperTreeNode<3>* node = tree->GetNode( index );
          index = node->GetChild( child );
          isLeaf = node->IsChildLeaf( child );
          break;
          }
        default:
          vtkGenericWarningMacro( "Bad branching factor " << this->GetBranchFactor() );
          return;
        } // case 1
      break;
    }
}
