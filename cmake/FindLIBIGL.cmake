# - Try to find the LIBIGL library
# Once done this will define
#
#  LIBIGL_FOUND - system has LIBIGL
#  LIBIGL_INCLUDE_DIR - **the** LIBIGL include directory
if(LIBIGL_FOUND)
    return()
endif()
find_path(LIBIGL_INCLUDE_DIR igl/readOBJ.h
    HINTS
        ${LIBIGL_DIR}
        ENV LIBIGL_DIR
    PATHS
        ${PROJECT_SOURCE_DIR}/external/libigl
        ${PROJECT_SOURCE_DIR}/../external/libigl
        ${PROJECT_SOURCE_DIR}/../../external/libigl
    PATH_SUFFIXES include
    NO_DEFAULT_PATH
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(LIBIGL
    "\nlibigl not found --- You can download it using:\n\tgit clone https://github.com/libigl/libigl.git ${CMAKE_SOURCE_DIR}/../libigl"
    LIBIGL_INCLUDE_DIR)
mark_as_advanced(LIBIGL_INCLUDE_DIR)

set(LIBIGL_WITH_OPENGL OFF CACHE INTERNAL "turn off OPENGL in LIBIGL")
set(LIBIGL_WITH_OPENGL_GLFW OFF CACHE INTERNAL "turn off OPENGL GLFW in LIBIGL")
set(LIBIGL_USE_STATIC_LIBRARY OFF CACHE INTERNAL "prefer STATIC build")

list(APPEND CMAKE_MODULE_PATH "${LIBIGL_INCLUDE_DIR}/../cmake")
message(STATUS "USE DIR: ${LIBIGL_INCLUDE_DIR}")
set(EIGEN_INCLUDE_DIR "${LIBIGL_INCLUDE_DIR}/../external/eigen")
include(libigl)