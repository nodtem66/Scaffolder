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
        ${PROJECT_SOURCE_DIR}/external/tbb
        ${PROJECT_SOURCE_DIR}/../external/tbb
        ${PROJECT_SOURCE_DIR}/../../external/tbb
    PATH_SUFFIXES include
    NO_DEFAULT_PATH
)
set(TBB_ROOT_DIR ${TBB_INCLUDE_DIR}/../)

set(TBB_TEST OFF CACHE INTERNAL "turn off")
set(TBB_EXAMPLES OFF CACHE INTERNAL "turn off")
set(TBB_STRICT OFF CACHE INTERNAL "turn off")
set(BUILD_SHARED_LIBS off CACHE INTERNAL "turn off dynamic LINK_DEPENDS (not recommended by intel)")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(TBB
    "\nTBB not found --- You can download it using:\n\tgit clone https://github.com/TBB/TBB ${CMAKE_SOURCE_DIR}/../TBB"
    TBB_INCLUDE_DIR)
mark_as_advanced(TBB_INCLUDE_DIR)

list(APPEND CMAKE_MODULE_PATH "${TBB_ROOT_DIR}")
message(STATUS "USE DIR: ${TBB_ROOT_DIR}")
add_subdirectory(${TBB_ROOT_DIR})