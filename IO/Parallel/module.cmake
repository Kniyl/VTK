vtk_module(vtkIOParallel
  GROUPS
    MPI
  DEPENDS
    vtkParallelCore
    vtkFiltersParallel
    vtkIOParallelMPI
    vtkIONetCDF
    vtkexodusII
  TEST_DEPENDS
    vtkTestingCore
  )
