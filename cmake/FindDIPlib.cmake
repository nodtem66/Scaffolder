# - Try to find the DIPlib library
# Once done this will define
#
#  DIPlib_FOUND - system has DIPlib
#  DIPlib_INCLUDE_DIR - **the** DIPlib include directory
if(DIPlib_FOUND)
    return()
endif()

find_path(DIPlib_INCLUDE_DIR diplib.h
    HINTS
        ${DIPlib_DIR}
        ENV DIPlib_DIR
    PATHS
        ${CMAKE_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/..
        ${CMAKE_SOURCE_DIR}/../..
        ${CMAKE_SOURCE_DIR}/diplib
        ${CMAKE_SOURCE_DIR}/../diplib
        ${CMAKE_SOURCE_DIR}/../../diplib        
        /usr
        /usr/local
        /usr/local/diplib
    PATH_SUFFIXES include
)
# TODO: Rewrite this cmake to add the source without compiled lib
set(DIPlib_LIB_DIR ${DIPlib_INCLUDE_DIR}/../lib)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(DIPlib
    "\nDIPlib not found --- You can download it using:\n\tgit clone https://github.com/DIPlib/diplib ${CMAKE_SOURCE_DIR}/../diplib"
    DIPlib_INCLUDE_DIR)
mark_as_advanced(DIPlib_INCLUDE_DIR)

message(STATUS "USE DIR: ${DIPlib_INCLUDE_DIR}/")