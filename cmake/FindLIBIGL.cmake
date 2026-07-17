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
set(LIBIGL_SKIP_DOWNLOAD ON CACHE INTERNAL "Skip download libigl")
# set(HUNTER_ENABLED ON CACHE INTERNAL "Enable Hunter package manager support")

message(STATUS "USE IGL DIR: ${LIBIGL_INCLUDE_DIR}")
# Suppress Eigen install rules — we don't want Eigen headers/libs in our install tree.
set(SCAFFOLDER_SKIP_EIGEN_INSTALL ON)
add_subdirectory("${LIBIGL_INCLUDE_DIR}/../")
message(STATUS "Eigen DIR: ${EIGEN3_INCLUDE_DIR}")