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
        ${CMAKE_SOURCE_DIR}/external/libigl
    PATH_SUFFIXES include
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(LIBIGL
    "\nlibigl not found --- You can download it using:\n\tgit clone https://github.com/libigl/libigl.git ${CMAKE_SOURCE_DIR}/../libigl"
    LIBIGL_INCLUDE_DIR)
mark_as_advanced(LIBIGL_INCLUDE_DIR)

set(LIBIGL_WITH_OPENGL OFF CACHE INTERNAL "")
set(LIBIGL_WITH_OPENGL_GLFW OFF CACHE INTERNAL "")

list(APPEND CMAKE_MODULE_PATH "${LIBIGL_INCLUDE_DIR}/../cmake")
include(libigl)
if (VCG_INCLUDE_DIR)
    set(eigen_version "")
    if (EXISTS "${VCG_INCLUDE_DIR}/eigenlib/eigen_version.txt")
        file(READ "${VCG_INCLUDE_DIR}/eigenlib/eigen_version.txt" eigen_version)
    endif()
    if (NOT eigen_version STREQUAL LIBIGL_EIGEN_VERSION)
        file(COPY "${LIBIGL_EXTERNAL}/eigen/Eigen" "${LIBIGL_EXTERNAL}/eigen/unsupported" DESTINATION "${VCG_INCLUDE_DIR}/eigenlib/")
        file(WRITE "${VCG_INCLUDE_DIR}/eigenlib/eigen_version.txt" "${LIBIGL_EIGEN_VERSION}")
    endif()
endif()