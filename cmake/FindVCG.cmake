# - Try to find the VCG library
# Once done this will define
#
#  VCG_FOUND - system has VCG
#  VCG_INCLUDE_DIR - **the** VCG include directory
if(VCG_FOUND)
    return()
endif()

find_path(VCG_INCLUDE_DIR vcg/complex/base.h
    HINTS
        ${VCG_DIR}
        ENV VCG_DIR
    PATHS
        ${PROJECT_SOURCE_DIR}/external/vcglib
        ${PROJECT_SOURCE_DIR}/../external/vcglib
        ${PROJECT_SOURCE_DIR}/../../external/vcglib
    PATH_SUFFIXES include
    NO_DEFAULT_PATH
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(VCG
    "\nVCG not found --- You can download it using:\n\tgit clone https://github.com/cnr-isti-vclab/vcglib.git ${CMAKE_SOURCE_DIR}/../vcglib"
    VCG_INCLUDE_DIR)
mark_as_advanced(VCG_INCLUDE_DIR)

message(STATUS "USE DIR: ${VCG_INCLUDE_DIR}/")
if (VCG_INCLUDE_DIR AND LIBIGL_INCLUDE_DIR)
    set(eigen_version "")
    if (EXISTS "${VCG_INCLUDE_DIR}/eigenlib/eigen_version.txt")
        file(READ "${VCG_INCLUDE_DIR}/eigenlib/eigen_version.txt" eigen_version)
    endif()
    if (NOT eigen_version STREQUAL LIBIGL_EIGEN_VERSION)
        file(COPY "${LIBIGL_EXTERNAL}/eigen/Eigen" "${LIBIGL_EXTERNAL}/eigen/unsupported" DESTINATION "${VCG_INCLUDE_DIR}/eigenlib/")
        file(WRITE "${VCG_INCLUDE_DIR}/eigenlib/eigen_version.txt" "${LIBIGL_EIGEN_VERSION}")
    endif()
endif()