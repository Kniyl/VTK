# Generate the VTKConfig.cmake file in the build tree.  Also configure
# one for installation.  The file tells external projects how to use
# VTK.

#-----------------------------------------------------------------------------
# Settings shared between the build tree and install tree.

IF(VTK_USE_MPI)
  SET(VTK_MPIRUN_EXE_CONFIG ${VTK_MPIRUN_EXE})
  SET(VTK_MPI_MAX_NUMPROCS_CONFIG ${VTK_MPI_MAX_NUMPROCS})
  SET(VTK_MPI_POSTFLAGS_CONFIG ${VTK_MPI_POSTFLAGS})
  SET(VTK_MPI_PREFLAGS_CONFIG ${VTK_MPI_PREFLAGS})
  SET(VTK_MPI_CLIENT_PREFLAGS_CONFIG ${VTK_MPI_CLIENT_PREFLAGS})
  SET(VTK_MPI_SERVER_PREFLAGS_CONFIG ${VTK_MPI_SERVER_PREFLAGS})
ELSE(VTK_USE_MPI)
  SET(VTK_MPIRUN_EXE_CONFIG "")
  SET(VTK_MPI_MAX_NUMPROCS_CONFIG "")
  SET(VTK_MPI_POSTFLAGS_CONFIG "")
  SET(VTK_MPI_PREFLAGS_CONFIG "")
  SET(VTK_MPI_CLIENT_PREFLAGS_CONFIG "")
  SET(VTK_MPI_SERVER_PREFLAGS_CONFIG "")
ENDIF(VTK_USE_MPI)

# Some setting for the remote testing
SET(VTK_MPI_CLIENT_POSTFLAGS_CONFIG ${VTK_MPI_CLIENT_POSTFLAGS})
SET(VTK_MPI_SERVER_POSTFLAGS_CONFIG ${VTK_MPI_SERVER_POSTFLAGS})

#-----------------------------------------------------------------------------
# Settings specific to the build tree.

# The "use" file.
SET(VTK_USE_FILE ${VTK_BINARY_DIR}/UseVTK.cmake)

# The build settings file.
SET(VTK_BUILD_SETTINGS_FILE ${VTK_BINARY_DIR}/VTKBuildSettings.cmake)

# The directory containing class list files for each kit.
SET(VTK_KITS_DIR_CONFIG ${VTK_BINARY_DIR}/Utilities)

# The wrapping hints file.
SET(VTK_WRAP_HINTS_CONFIG ${VTK_WRAP_HINTS})

# Library directory.
SET(VTK_LIBRARY_DIRS_CONFIG ${LIBRARY_OUTPUT_PATH})

# Runtime directory.
SET(VTK_RUNTIME_DIRS_CONFIG ${LIBRARY_OUTPUT_PATH})

# Determine the include directories needed.
SET(VTK_INCLUDE_DIRS_CONFIG
  ${VTK_INCLUDE_DIRS_BUILD_TREE}
  ${VTK_INCLUDE_DIRS_SOURCE_TREE}
  ${VTK_INCLUDE_DIRS_SYSTEM}
)

# Executable locations.
SET(VTK_TCL_HOME_CONFIG "")
SET(VTK_JAVA_JAR_CONFIG "")
SET(VTK_PARSE_JAVA_EXE_CONFIG "")
SET(VTK_WRAP_JAVA_EXE_CONFIG "")
SET(VTK_WRAP_PYTHON_EXE_CONFIG "")
SET(VTK_WRAP_PYTHON_INIT_EXE_CONFIG "")
SET(VTK_WRAP_TCL_EXE_CONFIG "")
SET(VTK_WRAP_TCL_INIT_EXE_CONFIG "")
IF(VTK_WRAP_TCL)
  SET(VTK_WRAP_TCL_EXE_CONFIG ${VTK_WRAP_TCL_EXE})
  SET(VTK_WRAP_TCL_INIT_EXE_CONFIG ${VTK_WRAP_TCL_INIT_EXE})
  SET(VTK_TCL_HOME_CONFIG ${VTK_TCL_HOME})
ENDIF(VTK_WRAP_TCL)
IF(VTK_WRAP_PYTHON)
  SET(VTK_WRAP_PYTHON_EXE_CONFIG ${VTK_WRAP_PYTHON_EXE})
  SET(VTK_WRAP_PYTHON_INIT_EXE_CONFIG ${VTK_WRAP_PYTHON_INIT_EXE})
ENDIF(VTK_WRAP_PYTHON)
IF(VTK_WRAP_JAVA)
  SET(VTK_PARSE_JAVA_EXE_CONFIG ${VTK_PARSE_JAVA_EXE})
  SET(VTK_WRAP_JAVA_EXE_CONFIG ${VTK_WRAP_JAVA_EXE})
  SET(VTK_JAVA_JAR_CONFIG ${LIBRARY_OUTPUT_PATH}/vtk.jar)
ENDIF(VTK_WRAP_JAVA)

# VTK style script locations.
SET(VTK_DOXYGEN_HOME_CONFIG ${VTK_SOURCE_DIR}/Utilities/Doxygen)
SET(VTK_HEADER_TESTING_PY_CONFIG ${VTK_SOURCE_DIR}/Common/Testing/HeaderTesting.py)
SET(VTK_FIND_STRING_TCL_CONFIG ${VTK_SOURCE_DIR}/Common/Testing/Tcl/FindString.tcl)
SET(VTK_PRINT_SELF_CHECK_TCL_CONFIG ${VTK_SOURCE_DIR}/Common/Testing/Tcl/PrintSelfCheck.tcl)
SET(VTK_RT_IMAGE_TEST_TCL_CONFIG ${VTK_SOURCE_DIR}/Common/Testing/Tcl/rtImageTest.tcl)

IF(VTK_USE_PARALLEL)
  SET(VTK_PRT_IMAGE_TEST_TCL_CONFIG ${VTK_SOURCE_DIR}/Common/Testing/Tcl/prtImageTest.tcl)
ELSE(VTK_USE_PARALLEL)
  SET(VTK_PRT_IMAGE_TEST_TCL_CONFIG "")
ENDIF(VTK_USE_PARALLEL)

# Location of tk internal headers provided by VTK.
IF(VTK_RENDERING_NEED_TK_INTERNAL)
  SET(VTK_TK_INTERNAL_DIR_CONFIG ${TK_INTERNAL_PATH})
ELSE(VTK_RENDERING_NEED_TK_INTERNAL)
  SET(VTK_TK_INTERNAL_DIR_CONFIG "")
ENDIF(VTK_RENDERING_NEED_TK_INTERNAL)

# CMake extension module directory.
SET(VTK_LOAD_CMAKE_EXTENSIONS_MACRO_CONFIG
    "${VTK_SOURCE_DIR}/CMake/vtkLoadCMakeExtensions.cmake")
SET(VTK_CMAKE_DIR_CONFIG
    "${VTK_SOURCE_DIR}/CMake")
SET(VTK_TCL_TK_MACROS_MODULE_CONFIG
    "${VTK_SOURCE_DIR}/CMake/vtkTclTkMacros.cmake")
SET(VTK_CMAKE_EXTENSIONS_DIR_CONFIG ${VTK_BINARY_DIR}/CMake)

# Library dependencies file.
SET(VTK_LIBRARY_DEPENDS_FILE "${VTK_BINARY_DIR}/VTKLibraryDepends.cmake")

# Build configuration information.
SET(VTK_CONFIGURATION_TYPES_CONFIG ${CMAKE_CONFIGURATION_TYPES})
SET(VTK_BUILD_TYPE_CONFIG ${CMAKE_BUILD_TYPE})

# Hack to give source tree access for a build tree configuration.
STRING(ASCII 35 VTK_STRING_POUND)
STRING(ASCII 64 VTK_STRING_AT)
SET(VTK_CONFIG_BACKWARD_COMPATIBILITY_HACK
  "\n${VTK_STRING_POUND} For backward compatability.  DO NOT USE.\nSET(VTK_SOURCE_DIR \"${VTK_SOURCE_DIR}\")\nIF(NOT TCL_LIBRARY)\n  SET(TCL_LIBRARY \"${TCL_LIBRARY}\" CACHE FILEPATH \"Location of Tcl library imported from VTK.  This may mean your project is depending on VTK to get this setting.  Consider using FindTCL.cmake.\")\nENDIF(NOT TCL_LIBRARY)\nIF(NOT TK_LIBRARY)\n  SET(TK_LIBRARY \"${TK_LIBRARY}\" CACHE FILEPATH \"Location of Tk library imported from VTK.  This may mean your project is depending on VTK to get this setting.  Consider using FindTCL.cmake.\")\nENDIF(NOT TK_LIBRARY)\nMARK_AS_ADVANCED(TCL_LIBRARY TK_LIBRARY)\n")

#-----------------------------------------------------------------------------
# Configure VTKConfig.cmake for the build tree.
CONFIGURE_FILE(${VTK_SOURCE_DIR}/VTKConfig.cmake.in
               ${VTK_BINARY_DIR}/VTKConfig.cmake @ONLY IMMEDIATE)

#-----------------------------------------------------------------------------
# Settings specific to the install tree.

# The "use" file.
SET(VTK_USE_FILE ${CMAKE_INSTALL_PREFIX}${VTK_INSTALL_PACKAGE_DIR}/UseVTK.cmake)

# The build settings file.
SET(VTK_BUILD_SETTINGS_FILE ${CMAKE_INSTALL_PREFIX}${VTK_INSTALL_PACKAGE_DIR}/VTKBuildSettings.cmake)

# The directory containing class list files for each kit.
SET(VTK_KITS_DIR_CONFIG ${CMAKE_INSTALL_PREFIX}${VTK_INSTALL_PACKAGE_DIR})

# The wrapping hints file.
IF(VTK_WRAP_HINTS)
  GET_FILENAME_COMPONENT(VTK_HINTS_FNAME ${VTK_WRAP_HINTS} NAME)
  SET(VTK_WRAP_HINTS_CONFIG ${CMAKE_INSTALL_PREFIX}${VTK_INSTALL_PACKAGE_DIR}/${VTK_HINTS_FNAME})
ENDIF(VTK_WRAP_HINTS)

# Include directories.
SET(VTK_INCLUDE_DIRS_CONFIG
  ${CMAKE_INSTALL_PREFIX}${VTK_INSTALL_INCLUDE_DIR}
  ${VTK_INCLUDE_DIRS_SYSTEM}
)

# Link directories.
SET(VTK_LIBRARY_DIRS_CONFIG ${CMAKE_INSTALL_PREFIX}${VTK_INSTALL_LIB_DIR})

# Runtime directories.
IF(WIN32)
  SET(VTK_RUNTIME_DIRS_CONFIG ${CMAKE_INSTALL_PREFIX}${VTK_INSTALL_BIN_DIR})
ELSE(WIN32)
  SET(VTK_RUNTIME_DIRS_CONFIG ${CMAKE_INSTALL_PREFIX}${VTK_INSTALL_LIB_DIR})
ENDIF(WIN32)

# Executable locations.
SET(VTK_TCL_HOME_CONFIG "")
SET(VTK_JAVA_JAR_CONFIG "")
SET(VTK_PARSE_JAVA_EXE_CONFIG "")
SET(VTK_WRAP_JAVA_EXE_CONFIG "")
SET(VTK_WRAP_PYTHON_EXE_CONFIG "")
SET(VTK_WRAP_PYTHON_INIT_EXE_CONFIG "")
SET(VTK_WRAP_TCL_EXE_CONFIG "")
SET(VTK_WRAP_TCL_INIT_EXE_CONFIG "")
SET(VTK_DOXYGEN_HOME_CONFIG "")
IF(VTK_WRAP_TCL)
  SET(VTK_WRAP_TCL_EXE_CONFIG ${CMAKE_INSTALL_PREFIX}${VTK_INSTALL_BIN_DIR}/vtkWrapTcl)
  SET(VTK_WRAP_TCL_INIT_EXE_CONFIG ${CMAKE_INSTALL_PREFIX}${VTK_INSTALL_BIN_DIR}/vtkWrapTclInit)
  SET(VTK_TCL_HOME_CONFIG ${CMAKE_INSTALL_PREFIX}${VTK_INSTALL_TCL_DIR}/tcl)
ENDIF(VTK_WRAP_TCL)
IF(VTK_WRAP_PYTHON)
  SET(VTK_WRAP_PYTHON_EXE_CONFIG ${CMAKE_INSTALL_PREFIX}${VTK_INSTALL_BIN_DIR}/vtkWrapPython)
  SET(VTK_WRAP_PYTHON_INIT_EXE_CONFIG ${CMAKE_INSTALL_PREFIX}${VTK_INSTALL_BIN_DIR}/vtkWrapPythonInit)
ENDIF(VTK_WRAP_PYTHON)
IF(VTK_WRAP_JAVA)
  SET(VTK_PARSE_JAVA_EXE_CONFIG ${CMAKE_INSTALL_PREFIX}${VTK_INSTALL_BIN_DIR}/vtkParseJava)
  SET(VTK_WRAP_JAVA_EXE_CONFIG ${CMAKE_INSTALL_PREFIX}${VTK_INSTALL_BIN_DIR}/vtkWrapJava)
  SET(VTK_JAVA_JAR_CONFIG ${CMAKE_INSTALL_PREFIX}${VTK_INSTALL_JAVA_DIR}/vtk.jar)
ENDIF(VTK_WRAP_JAVA)

# VTK style script locations.
SET(VTK_DOXYGEN_HOME_CONFIG ${CMAKE_INSTALL_PREFIX}${VTK_INSTALL_DOXYGEN_DIR})
SET(VTK_HEADER_TESTING_PY_CONFIG ${CMAKE_INSTALL_PREFIX}${VTK_INSTALL_PACKAGE_DIR}/testing/HeaderTesting.py)
SET(VTK_FIND_STRING_TCL_CONFIG ${CMAKE_INSTALL_PREFIX}${VTK_INSTALL_PACKAGE_DIR}/testing/FindString.tcl)
SET(VTK_PRINT_SELF_CHECK_TCL_CONFIG ${CMAKE_INSTALL_PREFIX}${VTK_INSTALL_PACKAGE_DIR}/testing/PrintSelfCheck.tcl)
SET(VTK_RT_IMAGE_TEST_TCL_CONFIG ${CMAKE_INSTALL_PREFIX}${VTK_INSTALL_PACKAGE_DIR}/testing/rtImageTest.tcl)

IF(VTK_USE_PARALLEL)
  SET(VTK_PRT_IMAGE_TEST_TCL_CONFIG ${CMAKE_INSTALL_PREFIX}${VTK_INSTALL_PACKAGE_DIR}/testing/prtImageTest.tcl)
ELSE(VTK_USE_PARALLEL)
  SET(VTK_PRT_IMAGE_TEST_TCL_CONFIG "")
ENDIF(VTK_USE_PARALLEL)

# Location of tk internal headers provided by VTK.
SET(VTK_TK_INTERNAL_DIR_CONFIG "")
IF(VTK_RENDERING_NEED_TK_INTERNAL AND TK_INTERNAL_PATH)
  VTK_GET_TCL_TK_VERSION ("TCL_TK_MAJOR_VERSION" "TCL_TK_MINOR_VERSION")
  SET (TCL_TK_VERSION "${TCL_TK_MAJOR_VERSION}.${TCL_TK_MINOR_VERSION}")
  IF("${TK_INTERNAL_PATH}" MATCHES 
     "Utilities/TclTk/internals/tk${TCL_TK_VERSION}")
    SET(VTK_TK_INTERNAL_DIR_CONFIG 
        ${CMAKE_INSTALL_PREFIX}${VTK_INSTALL_INCLUDE_DIR}/TclTk/internals/${TCL_TK_VERSION})
  ENDIF("${TK_INTERNAL_PATH}" MATCHES 
        "Utilities/TclTk/internals/tk${TCL_TK_VERSION}")
ENDIF(VTK_RENDERING_NEED_TK_INTERNAL AND TK_INTERNAL_PATH)

SET(VTK_TK_INTERNAL_ROOT_CONFIG ${CMAKE_INSTALL_PREFIX}${VTK_INSTALL_INCLUDE_DIR})

# CMake extension module directory and macro file.
SET(VTK_LOAD_CMAKE_EXTENSIONS_MACRO_CONFIG
    "${CMAKE_INSTALL_PREFIX}${VTK_INSTALL_PACKAGE_DIR}/CMake/vtkLoadCMakeExtensions.cmake")
SET(VTK_CMAKE_DIR_CONFIG
    "${CMAKE_INSTALL_PREFIX}${VTK_INSTALL_PACKAGE_DIR}/CMake")
SET(VTK_TCL_TK_MACROS_MODULE_CONFIG
    "${CMAKE_INSTALL_PREFIX}${VTK_INSTALL_PACKAGE_DIR}/CMake/vtkTclTkMacros.cmake")
SET(VTK_CMAKE_EXTENSIONS_DIR_CONFIG ${CMAKE_INSTALL_PREFIX}${VTK_INSTALL_PACKAGE_DIR}/CMake)

# Library dependencies file.
SET(VTK_LIBRARY_DEPENDS_FILE "${CMAKE_INSTALL_PREFIX}${VTK_INSTALL_PACKAGE_DIR}/VTKLibraryDepends.cmake")

# No backward compatibility hack needed for installed path
SET(VTK_CONFIG_BACKWARD_COMPATIBILITY_HACK)

#-----------------------------------------------------------------------------
# Configure VTKConfig.cmake for the install tree.
SET(VTK_CONFIGURATION_TYPES_CONFIG)
IF(CMAKE_CONFIGURATION_TYPES)
  # There are multiple build configurations.  Configure one
  # VTKConfig.cmake for each configuration.
  FOREACH(config ${CMAKE_CONFIGURATION_TYPES})
    SET(VTK_BUILD_TYPE_CONFIG ${config})
    CONFIGURE_FILE(${VTK_SOURCE_DIR}/VTKConfig.cmake.in
                   ${VTK_BINARY_DIR}/Utilities/${config}/VTKConfig.cmake
                   @ONLY IMMEDIATE)
  ENDFOREACH(config)

  # Install the config file corresponding to the build configuration
  # specified when building the install target.  The BUILD_TYPE variable
  # will be set while CMake is processing the install files.
  SET(DOLLAR "$")
  IF(NOT VTK_INSTALL_NO_DEVELOPMENT)
    INSTALL_FILES(${VTK_INSTALL_PACKAGE_DIR} FILES
      ${VTK_BINARY_DIR}/Utilities/${DOLLAR}{BUILD_TYPE}/VTKConfig.cmake)
  ENDIF(NOT VTK_INSTALL_NO_DEVELOPMENT)
ELSE(CMAKE_CONFIGURATION_TYPES)
  # There is only one build configuration.  Configure one VTKConfig.cmake.
  SET(VTK_BUILD_TYPE_CONFIG ${CMAKE_BUILD_TYPE})
  CONFIGURE_FILE(${VTK_SOURCE_DIR}/VTKConfig.cmake.in
                 ${VTK_BINARY_DIR}/Utilities/VTKConfig.cmake @ONLY IMMEDIATE)

  # Setup an install rule for the config file.
  IF(NOT VTK_INSTALL_NO_DEVELOPMENT)
    INSTALL_FILES(${VTK_INSTALL_PACKAGE_DIR} FILES
      ${VTK_BINARY_DIR}/Utilities/VTKConfig.cmake)
  ENDIF(NOT VTK_INSTALL_NO_DEVELOPMENT)
ENDIF(CMAKE_CONFIGURATION_TYPES)
