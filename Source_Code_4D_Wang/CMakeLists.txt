project(XTEngine)


set(PROJECT_PREFIX "XTEngine" CACHE STRING "The Prefix to be used for Preprocessor definitions")
set(PROJECT_INSTALL_HEADERS 0)
#-------------------------------------------------------------------------------
# Include the cmp project to perform all the basic configuration tests for the
# project.
#-------------------------------------------------------------------------------
set(CMP_HEADER_DIR ${XTEngine_BINARY_DIR})
set(CMP_CONFIGURATION_FILE_NAME "XTEngineConfiguration.h")
set(CMP_TYPES_FILE_NAME "XTEngineTypes.h")
set(CMP_VERSION_HEADER_FILE_NAME "XTEngineVersion.h")
set(CMP_VERSION_SOURCE_FILE_NAME "XTEngineVersion.cpp")
set(CMP_TOP_HEADER_FILE "XTEngine.h")
set(CMP_PROJECT_NAMESPACE "${PROJECT_NAME}")
set(CMP_PROJECT_NAME "${PROJECT_NAME}")
set(CMP_ENABLE_PLUGINS "0")
set(CMP_LIB_SEARCH_DIRS "")
set(CMP_GENERATE_VERSION_STRING 1)


INCLUDE (${CMP_SOURCE_DIR}/cmpProject.cmake)

cmpConfigureFileWithMD5Check(CONFIGURED_TEMPLATE_PATH ${XTEngine_SOURCE_DIR}/XTEngineConfiguration.h.in
                            GENERATED_FILE_PATH ${MuriReconstructionProj_BINARY_DIR}/XTEngine/${CMP_TOP_HEADER_FILE})

cmp_IDE_GENERATED_PROPERTIES("Generated"
                             "${CMP_HEADER_DIR}/${CMP_TOP_HEADER_FILE};${CMP_HEADER_DIR}/${CMP_TYPES_FILE_NAME};${CMP_HEADER_DIR}/${CMP_CONFIGURATION_FILE_NAME};${CMP_HEADER_DIR}/${CMP_VERSION_HEADER_FILE_NAME}"
                             "${CMP_HEADER_DIR}/${CMP_VERSION_SOURCE_FILE_NAME}")

#-------------------------------------------------------------------------------
# Set all the headers needed for this executable
#-------------------------------------------------------------------------------
set(XTEngine_HDRS
  ${XTEngine_SOURCE_DIR}/XT_AMatrix.h
  ${XTEngine_SOURCE_DIR}/XT_Constants.h
  ${XTEngine_SOURCE_DIR}/XT_Engine.h
  ${XTEngine_SOURCE_DIR}/XT_ICD_update.h
  ${XTEngine_SOURCE_DIR}/XT_IOMisc.h
  ${XTEngine_SOURCE_DIR}/XT_ImageProc.h
  ${XTEngine_SOURCE_DIR}/XT_Init.h
  ${XTEngine_SOURCE_DIR}/XT_MPI.h
  ${XTEngine_SOURCE_DIR}/XT_NHICD.h
  ${XTEngine_SOURCE_DIR}/XT_Profile.h
  ${XTEngine_SOURCE_DIR}/XT_genSinogram.h
  ${XTEngine_SOURCE_DIR}/allocate.h
  ${XTEngine_SOURCE_DIR}/randlib.h
  ${XTEngine_SOURCE_DIR}/tiff.h
)

#-------------------------------------------------------------------------------
# Set all the sources needed for this executable
#-------------------------------------------------------------------------------
set(XTEngine_SRCS
  ${XTEngine_SOURCE_DIR}/XT_AMatrix.c
  ${XTEngine_SOURCE_DIR}/XT_Engine.c
  ${XTEngine_SOURCE_DIR}/XT_ICD_update.c
  ${XTEngine_SOURCE_DIR}/XT_IOMisc.c
  ${XTEngine_SOURCE_DIR}/XT_ImageProc.c
  ${XTEngine_SOURCE_DIR}/XT_Init.c
  ${XTEngine_SOURCE_DIR}/XT_MPI.c
  ${XTEngine_SOURCE_DIR}/XT_NHICD.c
  ${XTEngine_SOURCE_DIR}/XT_Profile.c
  ${XTEngine_SOURCE_DIR}/XT_genSinogram.c
  ${XTEngine_SOURCE_DIR}/allocate.c
  ${XTEngine_SOURCE_DIR}/randlib.c
  ${XTEngine_SOURCE_DIR}/tiff.c
)
cmp_IDE_GENERATED_PROPERTIES("XTEngine" "${XTEngine_HDRS}" "${XTEngine_SRCS}")


#-------------------------------------------------------------------------------
# Group all the source and header files into a single CMake Variable
#-------------------------------------------------------------------------------
set(XTEngine_SOURCE_FILES
  ${CMP_HEADER_DIR}/${CMP_TOP_HEADER_FILE}
  ${CMP_HEADER_DIR}/${CMP_TYPES_FILE_NAME}
  ${CMP_HEADER_DIR}/${CMP_CONFIGURATION_FILE_NAME}
  ${CMP_HEADER_DIR}/${CMP_VERSION_HEADER_FILE_NAME}
  ${CMP_HEADER_DIR}/${CMP_VERSION_SOURCE_FILE_NAME}
  ${XTEngine_HDRS}
  ${XTEngine_SRCS}
)


#-------------------------------------------------------------------------------
# Add definitions for the executable
#-------------------------------------------------------------------------------
ADD_DEFINITIONS(-DPOSITIVITY_CONSTRAINT)
ADD_DEFINITIONS(-DNO_COST_CALCULATE)


#-------------------------------------------------------------------------------
# Add in our Include Directories
#-------------------------------------------------------------------------------
include_directories(${XTEngine_BINARY_DIR})

set(MATH_LIB "")
if(NOT MSVC)
  set(MATH_LIB "m")
endif()

# --------------------------------------------------------------------
# Setup the install rules for the various platforms
set(install_dir "tools")
if(WIN32)
    set(install_dir ".")
endif()


BuildToolBundle(
    TARGET XTEngine
    DEBUG_EXTENSION ${EXE_DEBUG_EXTENSION}
    VERSION_MAJOR ${DMuriReconstructionProj_VER_MAJOR}
    VERSION_MINOR ${MuriReconstructionProj_VER_MINOR}
    VERSION_PATCH ${MuriReconstructionProj_VER_PATCH}
    BINARY_DIR    ${${PROJECT_NAME}_BINARY_DIR}
    COMPONENT     Applications
    INSTALL_DEST  "${install_dir}"
    SOLUTION_FOLDER "Tools"
    SOURCES ${XTEngine_SOURCE_FILES}
    LINK_LIBRARIES  ${MATH_LIB} ${TIFF_LIBRARY} ${MPI_C_LIBRARIES}
)


