if(OpenCTM_FOUND)
    return()
endif()
find_path(OpenCTM_INCLUDE_DIR openctm.h
    HINTS
        ${OpenCTM_DIR}
        ENV OpenCTM_DIR
    PATHS
        ${PROJECT_SOURCE_DIR}/external/OpenCTM-1.0.3
        ${PROJECT_SOURCE_DIR}/../external/OpenCTM-1.0.3
        ${PROJECT_SOURCE_DIR}/../../external/OpenCTM-1.0.3
    PATH_SUFFIXES lib
    NO_DEFAULT_PATH
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(OpenCTM
    "\nOpenCTM not found --- You can download it and place into external/OpenCTM"
    OpenCTM_INCLUDE_DIR)
mark_as_advanced(OpenCTM_INCLUDE_DIR)

message(STATUS "USE DIR: ${OpenCTM_INCLUDE_DIR}")