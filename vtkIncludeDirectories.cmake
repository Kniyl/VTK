#-----------------------------------------------------------------------------
# Include directories for other projects installed on the system.
SET(VTK_INCLUDE_DIRS_SYSTEM "")
IF(VTK_USE_RENDERING)
  # OpenGL include directories.
  IF(APPLE)
    IF(VTK_USE_X)
       SET(VTK_INCLUDE_DIRS_SYSTEM ${VTK_INCLUDE_DIRS_SYSTEM}
           ${OPENGL_INCLUDE_DIR})
    ENDIF(VTK_USE_X)
  ELSE(APPLE)
       SET(VTK_INCLUDE_DIRS_SYSTEM ${VTK_INCLUDE_DIRS_SYSTEM}
           ${OPENGL_INCLUDE_DIR})
  ENDIF(APPLE)

  IF(VTK_USE_X)
    # X include directories.
    SET(VTK_INCLUDE_DIRS_SYSTEM ${VTK_INCLUDE_DIRS_SYSTEM}
        ${CMAKE_Xlib_INCLUDE_PATH} ${CMAKE_Xutil_INCLUDE_PATH})
  ENDIF(VTK_USE_X)

  IF(VTK_HAVE_VG500)
    # VolumePro VG 500 include directory.
    SET(VTK_INCLUDE_DIRS_SYSTEM ${VTK_INCLUDE_DIRS_SYSTEM}
        ${VLI_INCLUDE_PATH_FOR_VG500})
  ENDIF(VTK_HAVE_VG500)

  IF(VTK_HAVE_VP1000)
    # VolumePro VP 1000 include directory.
    SET(VTK_INCLUDE_DIRS_SYSTEM ${VTK_INCLUDE_DIRS_SYSTEM}
        ${VLI_INCLUDE_PATH_FOR_VP1000})
  ENDIF(VTK_HAVE_VP1000)

  IF(VTK_USE_MANGLED_MESA)
    # Mangled Mesa include directory.
    SET(VTK_INCLUDE_DIRS_SYSTEM ${VTK_INCLUDE_DIRS_SYSTEM}
        ${MESA_INCLUDE_PATH})
  ELSE(VTK_USE_MANGLED_MESA)
    # Off-Screen Mesa include directory.
    IF(VTK_OPENGL_HAS_OSMESA)
      IF(OSMESA_INCLUDE_DIR)
        SET(VTK_INCLUDE_DIRS_SYSTEM ${VTK_INCLUDE_DIRS_SYSTEM}
            ${OSMESA_INCLUDE_DIR})
      ENDIF(OSMESA_INCLUDE_DIR)
    ENDIF(VTK_OPENGL_HAS_OSMESA)
  ENDIF(VTK_USE_MANGLED_MESA)

  IF(VTK_USE_GL2PS)
    # gl2ps include directory.
    SET(VTK_INCLUDE_DIRS_SYSTEM ${VTK_INCLUDE_DIRS_SYSTEM} ${GL2PS_INCLUDE_PATH})
  ENDIF(VTK_USE_GL2PS)  

ENDIF(VTK_USE_RENDERING)

IF(VTK_USE_PARALLEL)
  IF(VTK_USE_MPI)
    # MPI include directory.
    SET(VTK_INCLUDE_DIRS_SYSTEM ${VTK_INCLUDE_DIRS_SYSTEM} ${MPI_INCLUDE_PATH})
  ENDIF(VTK_USE_MPI)
ENDIF(VTK_USE_PARALLEL)

IF(VTK_WRAP_TCL)
  SET(VTK_INCLUDE_DIRS_SYSTEM ${VTK_INCLUDE_DIRS_SYSTEM} ${TCL_INCLUDE_PATH})
ENDIF(VTK_WRAP_TCL)

IF(VTK_WRAP_PYTHON)
  # Python include directory.
  SET(VTK_INCLUDE_DIRS_SYSTEM ${VTK_INCLUDE_DIRS_SYSTEM}
      ${PYTHON_INCLUDE_PATH})
ENDIF(VTK_WRAP_PYTHON)

# VTK_INCLUDE_NEED_TK is set in toplevel CMakeLists.txt file.
IF(VTK_INCLUDE_NEED_TK)
  # Tcl/Tk include directories.
  SET(VTK_INCLUDE_DIRS_SYSTEM ${VTK_INCLUDE_DIRS_SYSTEM} ${TK_INCLUDE_PATH})
  IF (WIN32)
    SET(VTK_INCLUDE_DIRS_SYSTEM ${VTK_INCLUDE_DIRS_SYSTEM} ${TK_XLIB_PATH})
  ENDIF (WIN32)
ENDIF(VTK_INCLUDE_NEED_TK)

IF(VTK_WRAP_JAVA)
  # Java include directories.
  SET(VTK_INCLUDE_DIRS_SYSTEM ${VTK_INCLUDE_DIRS_SYSTEM}
      ${JAVA_INCLUDE_PATH} ${JAVA_INCLUDE_PATH2} ${JAVA_AWT_INCLUDE_PATH})
ENDIF(VTK_WRAP_JAVA)

#-----------------------------------------------------------------------------
# Include directories from the source tree.
SET(VTK_INCLUDE_DIRS_SOURCE_TREE "")
IF(VTK_USE_PARALLEL)
  SET(VTK_INCLUDE_DIRS_SOURCE_TREE ${VTK_INCLUDE_DIRS_SOURCE_TREE} ${VTK_SOURCE_DIR}/Parallel)
ENDIF(VTK_USE_PARALLEL)
IF(VTK_USE_HYBRID)
  SET(VTK_INCLUDE_DIRS_SOURCE_TREE ${VTK_INCLUDE_DIRS_SOURCE_TREE} ${VTK_SOURCE_DIR}/Hybrid)
ENDIF(VTK_USE_HYBRID)
IF(VTK_USE_PATENTED)
  SET(VTK_INCLUDE_DIRS_SOURCE_TREE ${VTK_INCLUDE_DIRS_SOURCE_TREE} ${VTK_SOURCE_DIR}/Patented)
ENDIF(VTK_USE_PATENTED)
IF(VTK_USE_RENDERING)
  SET(VTK_INCLUDE_DIRS_SOURCE_TREE ${VTK_INCLUDE_DIRS_SOURCE_TREE} ${VTK_SOURCE_DIR}/Rendering)
ENDIF(VTK_USE_RENDERING)

# These directories are always needed.
SET(VTK_INCLUDE_DIRS_SOURCE_TREE ${VTK_INCLUDE_DIRS_SOURCE_TREE}
  ${VTK_SOURCE_DIR}/IO
  ${VTK_SOURCE_DIR}/Imaging
  ${VTK_SOURCE_DIR}/Graphics
  ${VTK_SOURCE_DIR}/Filtering
  ${VTK_SOURCE_DIR}/Common
)



# Access to vtkRegressionTestImage.h.
SET(VTK_INCLUDE_DIRS_SOURCE_TREE ${VTK_INCLUDE_DIRS_SOURCE_TREE}
  ${VTK_SOURCE_DIR}/Common/Testing/Cxx
)

#-----------------------------------------------------------------------------
# Include directories from the build tree.
SET(VTK_INCLUDE_DIRS_BUILD_TREE ${VTK_BINARY_DIR})

#-----------------------------------------------------------------------------
# Include directories needed for .cxx files in VTK.  These include
# directories will NOT be available to user projects.
SET(VTK_INCLUDE_DIRS_BUILD_TREE_CXX
  ${VTKFREETYPE_SOURCE_DIR}/include
  ${VTKFREETYPE_BINARY_DIR}
  ${VTKFTGL_SOURCE_DIR}/src
  ${VTKFTGL_BINARY_DIR}
)

IF(VTK_RENDERING_NEED_TK_INTERNAL)
  # Need access to internal Tk headers for the vtkTk... widget .cxx files.
  SET(VTK_INCLUDE_DIRS_BUILD_TREE_CXX ${VTK_INCLUDE_DIRS_BUILD_TREE_CXX}
      ${TK_INTERNAL_PATH})
ENDIF(VTK_RENDERING_NEED_TK_INTERNAL)

IF (VTK_USE_MATROX_IMAGING)
  # Need access to mil.h include file for vtkMILVideoSource.cxx.
  SET(VTK_INCLUDE_DIRS_BUILD_TREE_CXX ${VTK_INCLUDE_DIRS_BUILD_TREE_CXX}
      ${MIL_INCLUDE_PATH})
ENDIF (VTK_USE_MATROX_IMAGING)

#-----------------------------------------------------------------------------
# Include directories for 3rd-party utilities provided by VTK.
VTK_THIRD_PARTY_INCLUDE(ZLIB zlib)
VTK_THIRD_PARTY_INCLUDE(JPEG jpeg)
VTK_THIRD_PARTY_INCLUDE(PNG  png)
VTK_THIRD_PARTY_INCLUDE(TIFF tiff)
VTK_THIRD_PARTY_INCLUDE(EXPAT expat)

# DICOM support is disabled until the following are fixed by the authors:
#  - Make it build with and without VTK_USE_ANSI_STDLIB.
#  - Remove all build warnings from dashboard.
#  - Add a test to bring up coverage of vtkDICOMImageReader.cxx
#    and actually read a DICOM file.
#IF(VTK_USE_ANSI_STDLIB)
#  VTK_THIRD_PARTY_INCLUDE(DICOMParser DICOMParser)
#ENDIF(VTK_USE_ANSI_STDLIB)
