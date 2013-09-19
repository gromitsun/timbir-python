project(MuriReconstructionProj)
cmake_minimum_required(VERSION 2.8.9)

include (CMakeParseArguments)

# ---------- Setup output Directories -------------------------
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY
  ${PROJECT_BINARY_DIR}/Bin
  CACHE PATH
  "Single Directory for all Libraries"
  )

# --------- Setup the Executable output Directory -------------
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY
  ${PROJECT_BINARY_DIR}/Bin
  CACHE PATH
  "Single Directory for all Executables."
  )

# --------- Setup the Executable output Directory -------------
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY
  ${PROJECT_BINARY_DIR}/Bin
  CACHE PATH
  "Single Directory for all static libraries."
  )


#Enable use of Solution Folders
set_property(GLOBAL PROPERTY USE_FOLDERS ON)


set(PROJECT_CODE_DIR ${PROJECT_SOURCE_DIR}/Source_Code_4D)
set(PROJECT_RESOURCES_DIR ${MuriReconstructionProj_SOURCE_DIR}/Resources)
set(EXE_DEBUG_EXTENSION "_debug")
set(PROJECT_INSTALL_HEADERS "0")
set(PROJECT_INSTALL_EXECUTABLES "1")
set(MuriReconstructionProj_SUPPORT_DIR ${PROJECT_SOURCE_DIR}/Support)

set(CMP_SOURCE_DIR ${MuriReconstructionProj_SOURCE_DIR}/Support/cmp)
INCLUDE (${CMP_SOURCE_DIR}/cmpCMakeMacros.cmake )


set(CMP_VERSION_MAJOR "1.0")


# --------------------------------------------------------------------
# Over ride CMake's built in module directory by prepending cmp's module
# directory first
set(CMAKE_MODULE_PATH ${CMP_MODULES_SOURCE_DIR} ${CMAKE_MODULE_PATH})

# Build shared libraries
OPTION (BUILD_SHARED_LIBS "Build Shared Libraries" OFF)
set(LIB_TYPE STATIC)
if(BUILD_SHARED_LIBS)
    set(LIB_TYPE SHARED)
    list(APPEND CMP_LIB_SEARCH_DIRS  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY} )
    if(MSVC)
        list(APPEND CMP_LIB_SEARCH_DIRS ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/Debug
                                        ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/Release )
    endif()
endif(BUILD_SHARED_LIBS)

# Include some directories for the compiler
include_directories(${PROJECT_BINARY_DIR})


# Prepend our own "FindTiff" code in front of what CMake itself provides
set (MuriReconstructionProj_TIFF_SUPPORT)
FIND_PACKAGE(Tiff)
IF (NOT TIFF_FOUND)
      MESSAGE (FATAL_ERROR " Tiff Library was NOT Found and is needed.")
ENDIF (NOT TIFF_FOUND)
set (MuriReconstructionProj_TIFF_SUPPORT "1")
INCLUDE_DIRECTORIES(${TIFF_INCLUDE_DIRS})
# Copy DLL libraries into the correct build location on MSVC platforms
if (HAVE_TIFF_DLL)
CMP_COPY_DEPENDENT_LIBRARIES ("tiff")
CMP_LIBRARIES_INSTALL_RULES ("tiff" "." )
endif()

# --------------------------------------------------------------------
# No HDF5 Support at this point
# --------------------------------------------------------------------
set(MuriReconstructionProj_HDF5_SUPPORT "0")

# --------------------------------------------------------------------
# Find the openmp Libraries as those are needed. Apple's compilers either
# do not support openmp (Clang) or generate broken code (GCC 4.2). Using
# Intel Threading Building Blocks (TBB) for this would be a better cross
# platform solution
# --------------------------------------------------------------------
set(MuriReconstructionProj_OpenMP_SUPPORT "0")
if(NOT APPLE)
	find_package(OpenMP)
	if (OPENMP_FOUND)
		set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
		set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
		set(MuriReconstructionProj_OpenMP_SUPPORT "1")
	endif()
endif()

# --------------------------------------------------------------------
# Find the MPI Libraries as those are needed
# --------------------------------------------------------------------
FIND_PACKAGE(MPI)
IF(NOT MPI_FOUND)
	message(FATAL_ERROR "MPI Libraries are needed for this project")
ENDIF()
INCLUDE_DIRECTORIES(${MPI_C_INCLUDE_PATH} ${MPI_CXX_INCLUDE_PATH})

# --------------------------------------------------------------------
# Turn off lots of MSVC Warnings
# --------------------------------------------------------------------
if(WIN32)
	option(MuriReconstructionProj_DISABLE_MSVC_WARNINGS "Disable some MSVC Compiler warnings" OFF)
	if(MuriReconstructionProj_DISABLE_MSVC_WARNINGS)
		add_definitions(-DMuriReconstruction_DISABLE_MSVC_WARNINGS)
	endif()
	add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif()


# --------------------------------------------------------------------
# Add in some compiler definitions
# --------------------------------------------------------------------
if( CMAKE_BUILD_TYPE MATCHES Debug )
  ADD_DEFINITIONS(-DDEBUG)
  if(NOT MSVC)
  ADD_DEFINITIONS(-Wall)
  endif()
ENDif( CMAKE_BUILD_TYPE MATCHES Debug )

# On linux we need to set this because some of the libraries are Static
# and some are shared.
if( CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" AND NOT MSVC )
  ADD_DEFINITIONS(-fPIC)
ENDif()


# --------------------------------------------------------------------
# If was are using GCC, make the compiler messages on a single line
if(CMAKE_COMPILER_IS_GNUCC)
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fmessage-length=0")
ENDif(CMAKE_COMPILER_IS_GNUCC)
if(CMAKE_COMPILER_IS_GNUCXX)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fmessage-length=0")
ENDif(CMAKE_COMPILER_IS_GNUCXX)


include_directories(${MuriReconstructionProj_SOURCE_DIR}/Source_Code_4D)
include_directories(${PROJECT_BINARY_DIR})

ADD_SUBDIRECTORY( ${PROJECT_SOURCE_DIR}/Source_Code_4D ${PROJECT_BINARY_DIR}/XTEngine)
