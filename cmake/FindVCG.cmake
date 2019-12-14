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
        ${CMAKE_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/..
        ${CMAKE_SOURCE_DIR}/../..
        ${CMAKE_SOURCE_DIR}/vcglib
        ${CMAKE_SOURCE_DIR}/../vcglib
        ${CMAKE_SOURCE_DIR}/../../vcglib        
        /usr
        /usr/local
        /usr/local/vcglib
    PATH_SUFFIXES include
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(VCG
    "\nVCG not found --- You can download it using:\n\tgit clone https://github.com/cnr-isti-vclab/vcglib.git ${CMAKE_SOURCE_DIR}/../vcglib"
    VCG_INCLUDE_DIR)
mark_as_advanced(VCG_INCLUDE_DIR)

message(STATUS "INCLUDE DIR: ${VCG_INCLUDE_DIR}/")
# include(VCG)