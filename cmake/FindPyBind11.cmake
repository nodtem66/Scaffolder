if(PyBind11_FOUND)
    return()
endif()
find_path(PyBind11_INCLUDE_DIR PyBind11/PyBind11.h
    HINTS
        ${PyBind11_DIR}
        ENV PyBind11_DIR
    PATHS
        ${PROJECT_SOURCE_DIR}/external/pybind11
        ${PROJECT_SOURCE_DIR}/../external/pybind11
        ${PROJECT_SOURCE_DIR}/../../external/pybind11
    PATH_SUFFIXES include
    NO_DEFAULT_PATH
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PyBind11
    "\nPyBind11 not found --- You can download it"
    PyBind11_INCLUDE_DIR)
mark_as_advanced(PyBind11_INCLUDE_DIR)

message(STATUS "USE DIR: ${PyBind11_INCLUDE_DIR}")
add_subdirectory("${PyBind11_INCLUDE_DIR}/../")