# - Try to find the TBB library
# Once done this will define
#
#  TBB_FOUND - system has TBB
#  TBB_INCLUDE_DIR - **the** TBB include directory
if(TBB_FOUND)
    return()
endif()

find_path(TBB_INCLUDE_DIR tbb/tbb.h
    HINTS
        ${TBB_DIR}
        ENV TBB_DIR
    PATHS
        ${CMAKE_SOURCE_DIR}/external/tbb    
    PATH_SUFFIXES include
)
set(TBB_ROOT_DIR ${TBB_INCLUDE_DIR}/../)

set(TBB_BUILD_SHARED OFF CACHE INTERNAL "")
set(TBB_BUILD_STATIC ON CACHE INTERNAL "")
set(TBB_BUILD_TESTS OFF CACHE INTERNAL "")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(TBB
    "\nTBB not found --- You can download it using:\n\tgit clone https://github.com/TBB/TBB ${CMAKE_SOURCE_DIR}/../TBB"
    TBB_INCLUDE_DIR)
mark_as_advanced(TBB_INCLUDE_DIR)

list(APPEND CMAKE_MODULE_PATH "${TBB_ROOT_DIR}")
message(STATUS "USE DIR: ${TBB_ROOT_DIR}")
add_subdirectory(${TBB_ROOT_DIR})