vtk_module(vtkIOParallel
  GROUPS
    StandAlone
  DEPENDS
    vtkParallelCore
    vtkFiltersParallel
    vtkIONetCDF
    vtkIOXML
    vtkIOImage
  PRIVATE_DEPENDS
    vtkexodusII
  TEST_DEPENDS
    vtkParallelMPI
    vtkRenderingParallel
    vtkTestingCore
    vtkTestingRendering
  )
