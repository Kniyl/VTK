/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkCommon.h
  Language:  C++
  Date:      $Date$
  Version:   $Revision$


Copyright (c) 1993-1999 Ken Martin, Will Schroeder, Bill Lorensen.

This software is copyrighted by Ken Martin, Will Schroeder and Bill Lorensen.
The following terms apply to all files associated with the software unless
explicitly disclaimed in individual files. This copyright specifically does
not apply to the related textbook "The Visualization Toolkit" ISBN
013199837-4 published by Prentice Hall which is covered by its own copyright.

The authors hereby grant permission to use, copy, and distribute this
software and its documentation for any purpose, provided that existing
copyright notices are retained in all copies and that this notice is included
verbatim in any distributions. Additionally, the authors grant permission to
modify this software and its documentation for any purpose, provided that
such modifications are not distributed without the explicit consent of the
authors and that existing copyright notices are retained in all copies. Some
of the algorithms implemented by this software are patented, observe all
applicable patent law.

IN NO EVENT SHALL THE AUTHORS OR DISTRIBUTORS BE LIABLE TO ANY PARTY FOR
DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES ARISING OUT
OF THE USE OF THIS SOFTWARE, ITS DOCUMENTATION, OR ANY DERIVATIVES THEREOF,
EVEN IF THE AUTHORS HAVE BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

THE AUTHORS AND DISTRIBUTORS SPECIFICALLY DISCLAIM ANY WARRANTIES, INCLUDING,
BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
PARTICULAR PURPOSE, AND NON-INFRINGEMENT.  THIS SOFTWARE IS PROVIDED ON AN
"AS IS" BASIS, AND THE AUTHORS AND DISTRIBUTORS HAVE NO OBLIGATION TO PROVIDE
MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.


=========================================================================*/

#include "vtkAbstractMapper.h"
#include "vtkActor2D.h"
#include "vtkActor2DCollection.h"
#include "vtkAttributeData.h"
#include "vtkBitArray.h"
#include "vtkByteSwap.h"
#include "vtkCell.h"
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkCellLinks.h"
#include "vtkCellTypes.h"
#include "vtkCharArray.h"
#include "vtkCollection.h"
#include "vtkContourValues.h"
#include "vtkCoordinate.h"
#include "vtkDataArray.h"
#include "vtkDataObject.h"
#include "vtkDataSet.h"
#include "vtkDataSetAttributes.h"
#include "vtkDataSetCollection.h"
#include "vtkDoubleArray.h"
#include "vtkEdgeTable.h"
#include "vtkEmptyCell.h"
#include "vtkFieldData.h"
#include "vtkFloatArray.h"
#include "vtkFloatNormals.h"
#include "vtkFloatPoints.h"
#include "vtkFloatScalars.h"
#include "vtkFloatTCoords.h"
#include "vtkFloatTensors.h"
#include "vtkFloatVectors.h"
#include "vtkGenericCell.h"
#include "vtkHexahedron.h"
#include "vtkIdList.h"
#include "vtkImageData.h"
#include "vtkImageSource.h"
#include "vtkImageToStructuredPoints.h"
#include "vtkImplicitFunction.h"
#include "vtkImplicitFunctionCollection.h"
#include "vtkIndent.h"
#include "vtkIntArray.h"
#include "vtkLine.h"
#include "vtkLocator.h"
#include "vtkLogLookupTable.h"
#include "vtkLongArray.h"
#include "vtkLookupTable.h"
#include "vtkMapper2D.h"
#include "vtkMath.h"
#include "vtkMatrix4x4.h"
#include "vtkMergePoints2D.h"
#include "vtkMultiThreader.h"
#include "vtkMutexLock.h"
#include "vtkNormals.h"
#include "vtkObject.h"
#include "vtkPixel.h"
#include "vtkPlane.h"
#include "vtkPointData.h"
#include "vtkPointLocator.h"
#include "vtkPointLocator2D.h"
#include "vtkPointSet.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkPolyLine.h"
#include "vtkPolyVertex.h"
#include "vtkPolygon.h"
#include "vtkPriorityQueue.h"
#include "vtkProcessObject.h"
#include "vtkProcessStatistics.h"
#include "vtkProp.h"
#include "vtkPropCollection.h"
#include "vtkProperty2D.h"
#include "vtkPyramid.h"
#include "vtkQuad.h"
#include "vtkRectilinearGrid.h"
#include "vtkObject.h"
#include "vtkScalars.h"
#include "vtkShortArray.h"
#include "vtkSource.h"
#include "vtkStack.h"
#include "vtkStructuredData.h"
#include "vtkStructuredGrid.h"
#include "vtkStructuredPoints.h"
#include "vtkTCoords.h"
#include "vtkTensors.h"
#include "vtkTetra.h"
#include "vtkTimeStamp.h"
#include "vtkTimerLog.h"
#include "vtkTransform.h"
#include "vtkTransformCollection.h"
#include "vtkTriangle.h"
#include "vtkTriangleStrip.h"
#include "vtkUnsignedCharArray.h"
#include "vtkUnsignedShortArray.h"
#include "vtkUnsignedIntArray.h"
#include "vtkUnsignedLongArray.h"
#include "vtkUnstructuredGrid.h"
#include "vtkVectors.h"
#include "vtkVersion.h"
#include "vtkVertex.h"
#include "vtkViewport.h"
#include "vtkVoidArray.h"
#include "vtkVoxel.h"
#include "vtkWedge.h"
#include "vtkWindow.h"
#include "vtkWindowLevelLookupTable.h"
#include "vtkWindowToImageFilter.h"
