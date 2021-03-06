set(VTK_SMP_IMPLEMENTATION_TYPE "Sequential" CACHE STRING
  "Name of the implementation of SMP techniques."
)
set_property(CACHE VTK_SMP_IMPLEMENTATION_TYPE PROPERTY STRINGS Sequential Kaapi OpenMP TBB)

if( NOT ("${VTK_SMP_IMPLEMENTATION_TYPE}" STREQUAL "Kaapi" OR
         "${VTK_SMP_IMPLEMENTATION_TYPE}" STREQUAL "TBB" OR
         "${VTK_SMP_IMPLEMENTATION_TYPE}" STREQUAL "OpenMP"))
  set(VTK_SMP_IMPLEMENTATION_TYPE "Sequential")
endif()

set(VTK_SMP_CMAKE_IMPLEMENTATION_FILE
  ${CMAKE_CURRENT_SOURCE_DIR}/${VTK_SMP_IMPLEMENTATION_TYPE}/vtkSMP${VTK_SMP_IMPLEMENTATION_TYPE}Implementation.cmake)
include(${VTK_SMP_CMAKE_IMPLEMENTATION_FILE})

message(STATUS "Module vtkParallelSMP uses the " ${VTK_SMP_IMPLEMENTATION_TYPE} " runtime")
set(vtkParallelSMP_SYSTEM_INCLUDE_DIRS ${VTK_SMP_IMPLEMENTATION_HEADER_DIR})

#internal implementation details
set(SMP_Internals
  vtkDummyMergeFunctor.cxx
  vtkLockPointMerger.cxx
  vtkOffsetManager.cxx
  vtkParallelCellMerger.cxx
  vtkParallelPointMerger.cxx
  vtkTask.cxx
  ${VTK_SMP_IMPLEMENTATION_TYPE}/vtkSMPImplementation.cxx
)

#components that developers use to make SMP filters from
set(SMP_Core
  vtkParallelTree.cxx #required by TraverseOp
  vtkLocalData.cxx
  vtkTreeFunctor.cxx #requirement of ForEachOp
  vtkTreeFunctorInitializable.cxx
  vtkRangeFunctor.cxx
  vtkRangeFunctorInitializable.cxx
  vtkRange.cxx
  vtkRange1D.cxx
  vtkParallelOperators.cxx
  vtkMergeDataSets.cxx
  vtkThreadLocal.cxx
)

#end user SMP accelerated vtk facilities (mostly filters)
set(SMP_Filters
  vtkParallelCompositeDataPipeline.cxx
  vtkSplittingAlgorithm.cxx
  vtkSMPMergePoints.cxx
  vtkSMPMinMaxTree.cxx
  vtkSMPTransform.cxx
  vtkSMPWarpVector.cxx
  vtkSMPZCurve.cxx
  vtkSMPClipDataSet.cxx
  vtkSMPContourFilter.cxx
)

set(Module_SRCS
  ${SMP_Internals}
  ${SMP_Core}
  ${SMP_Filters}
)

vtk_module_library(vtkParallelSMP ${Module_SRCS})
target_link_libraries(vtkParallelSMP ${VTK_SMP_IMPLEMENTATION_LIBRARIES})
