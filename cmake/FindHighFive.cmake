# - Try to find the HighFive library
# Once done this will define
#
#  HighFive_FOUND - system has HighFive
#  HighFive_INCLUDE_DIR - **the** HighFive include directory
if(HighFive_FOUND)
    return()
endif()
find_path(HighFive_INCLUDE_DIR H5Easy.hpp
    HINTS
        ${HighFive_DIR}
        ENV HighFive_DIR
    PATHS
        ${CMAKE_SOURCE_DIR}/HighFive
        ${CMAKE_SOURCE_DIR}/../HighFive
        ${CMAKE_SOURCE_DIR}/../../HighFive
        ${CMAKE_SOURCE_DIR}/../..
        ${CMAKE_SOURCE_DIR}/..
        ${CMAKE_SOURCE_DIR}
        /usr
        /usr/local
        /usr/local/igl/HighFive
    PATH_SUFFIXES include/highfive
)

# find_path(H5_INCLUDE_DIR H5Ppublic.h)
# LIST(APPEND HighFive_INCLUDE_DIR ${H5_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(HighFive
    "\nHighFive not found --- You can download it using:\n\tgit clone https://github.com/BlueBrain/HighFive ${CMAKE_SOURCE_DIR}/../HighFive"
    HighFive_INCLUDE_DIR)
mark_as_advanced(HighFive_INCLUDE_DIR)
message(STATUS "USE DIR: ${HighFive_INCLUDE_DIR}/")