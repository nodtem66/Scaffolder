cmake_minimum_required(VERSION 3.20)
cmake_policy(VERSION 3.20)
cmake_policy(SET CMP0135 NEW)
cmake_policy(SET CMP0074 NEW)
cmake_policy(SET CMP0169 OLD)

project(Scaffolder
  DESCRIPTION "Scaffolder is a tool for generating scaffolds for 3D printing."
  LANGUAGES C CXX
)

if(DEFINED VERSION)
  message(STATUS "Scaffolder version: ${VERSION}")
else()
  set(VERSION "dev")
  message(STATUS "Scaffolder version: ${VERSION} (development build)")
endif()

# -----------------------------------------------------------------------------
# Project layout and options
# -----------------------------------------------------------------------------
set(SCAFFOLDER_ROOT "${CMAKE_CURRENT_SOURCE_DIR}")
set(SCAFFOLDER_INCLUDE_DIR "${SCAFFOLDER_ROOT}/include")
set(SCAFFOLDER_SOURCE_DIR "${SCAFFOLDER_ROOT}/src")
set(SCAFFOLDER_CMAKE_DIR "${SCAFFOLDER_ROOT}/cmake")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_WARN_DEPRECATED OFF CACHE BOOL "" FORCE)

option(LIBIGL_WITH_OPENGL        "Use OpenGL"             OFF)
option(LIBIGL_WITH_OPENGL_GLFW   "Use GLFW"               OFF)
option(LIBIGL_BUILD_TESTS        "Build libigl unit test" OFF)
option(LIBIGL_BUILD_TUTORIALS    "Build libigl tutorial"  OFF)
option(LIBIGL_USE_STATIC_LIBRARY "Use LibIGL as static library" ON)
option(LIBIGL_SKIP_DOWNLOAD      "Skip download libigl"         ON)
option(IGL_STATIC_RUNTIME        "Use LibIGL as static library" ON)
option(SCAFFOLD_BUILD_BINARY     "Build binary"                 ON)
option(SCAFFOLD_BUILD_PYTHON     "Build python library"         OFF)
option(SCAFFOLDER_BUILD_WITH_OPENMP "USE OpenMP"                ON)

# -----------------------------------------------------------------------------
# External dependency
# -----------------------------------------------------------------------------
include("${SCAFFOLDER_CMAKE_DIR}/CPM.cmake")
set(CPM_SOURCE_CACHE "${SCAFFOLDER_ROOT}/.cpm_cache/")
CPMAddPackage("gh:jarro2783/cxxopts@3.3.1#44380e5a44706ab7347f400698c703eb2a196202")
CPMAddPackage(
  URI "https://www.lua.org/ftp/lua-5.4.8.tar.gz#SHA256=4f18ddae154e793e46eeab727c59ef1c0c0c2b744e7b94219710d76f530629ae"
  DOWNLOAD_ONLY YES
)
CPMAddPackage("gh:libigl/libigl@v2.6.0#40e7900ccbd767f1f360e0eb10f0f1a6432e0993")
CPMAddPackage(
  URI "gh:Danny02/OpenCTM#9021bc220f47693005699c1a69e7b404c0f39677"
  PATCHES "${SCAFFOLDER_CMAKE_DIR}/openctm_fix_min_version.patch"
)
CPMAddPackage("gh:pybind/pybind11@v2.13.6#a2e59f0e7065404b44dfe92a28aca47ba1378dc4")
CPMAddPackage(
  URI "gh:ThePhD/sol2@v3.3.0#eba86625b707e3c8c99bbfc4624e51f42dc9e561"
  PATCHES "${SCAFFOLDER_CMAKE_DIR}/sol2_fix_for_clang.patch"
)
CPMAddPackage("gh:cnr-isti-vclab/vcglib#c94ef4e12e9ea3ae986d9af91005be8328d13719") # tag 2025.07
CPMAddPackage("gh:cpm-cmake/CPMLicenses.cmake@0.0.7#ca42334d561b83e499b11cf55fe05d97a0767fb9")

# -----------------------------------------------------------------------------
# OpenMP configuration
# -----------------------------------------------------------------------------
if(SCAFFOLDER_BUILD_WITH_OPENMP)
  find_package(OpenMP REQUIRED)
  message(STATUS "OpenMP: ${OpenMP_CXX_FLAGS}")
endif()

# -----------------------------------------------------------------------------
# Common target setup
# -----------------------------------------------------------------------------
set(SCAFFOLDER_COMMON_INCLUDE_DIRS
  "${SCAFFOLDER_INCLUDE_DIR}"
  "${cxxopts_SOURCE_DIR}/include/"
  "${progress-cpp_SOURCE_DIR}/include/progresscpp/"
  "${lua_SOURCE_DIR}/src"
  "${vcglib_SOURCE_DIR}"
  "${libigl_SOURCE_DIR}/include/"
  "${OpenCTM_SOURCE_DIR}/lib/"
  "${OpenCTM_SOURCE_DIR}/lib/liblzma/"
  "${sol2_SOURCE_DIR}"
)

file(GLOB lua_sources_without_binary
  "${lua_SOURCE_DIR}/src/*.c"
)
list(REMOVE_ITEM lua_sources_without_binary "${lua_SOURCE_DIR}/src/lua.c" "${lua_SOURCE_DIR}/src/luac.c")
add_library(lua STATIC ${lua_sources_without_binary})
target_include_directories(lua PRIVATE "${lua_SOURCE_DIR}/src")

set(SCAFFOLDER_COMMON_LINK_LIBS Eigen3::Eigen lua sol2)

message(STATUS "SCAFFOLDER_COMMON_INCLUDE_DIRS = ${SCAFFOLDER_COMMON_INCLUDE_DIRS}")
message(STATUS "SCAFFOLDER_COMMON_LINK_LIBS = ${SCAFFOLDER_COMMON_LINK_LIBS}")

if(TARGET OpenMP::OpenMP_CXX)
  list(APPEND SCAFFOLDER_COMMON_LINK_LIBS OpenMP::OpenMP_CXX)
endif()

function(scaffolder_configure_target target_name)
  target_include_directories(${target_name} PRIVATE ${SCAFFOLDER_COMMON_INCLUDE_DIRS})
  target_link_libraries(${target_name} PRIVATE ${SCAFFOLDER_COMMON_LINK_LIBS} cxxopts)
  target_compile_options(${target_name} PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>: -Wall -Wextra -pedantic>
    $<$<CXX_COMPILER_ID:MSVC>: /MP /bigobj /D_CRT_SECURE_NO_WARNINGS /Zc:__cplusplus
    /wd4267 /wd4244 /wd4018 /wd4305 /wd4667 /wd4996 /wd4804>
  )
  target_compile_definitions(${target_name} PRIVATE
    OPENCTM_STATIC
    $<$<PLATFORM_ID:Apple>: _DARWIN_C_SOURCE _GNU_SOURCE _XOPEN_SOURCE=600>
    $<$<CXX_COMPILER_ID:MSVC>: NOMINMAX _USE_MATH_DEFINES WIN32_LEAN_AND_MEAN>
    $<IF:$<TARGET_EXISTS:OpenMP::OpenMP_CXX>,USE_PARALLEL=1,USE_PARALLEL=0>
  )
  if(DEFINED VERSION)
    target_compile_definitions(${target_name} PRIVATE SCAFFOLDER_VERSION="${VERSION}")
  endif()
endfunction()

function(write_third_party_license file_name)
  message(STATUS "CPM_PACKAGES = ${CPM_PACKAGES}")
  write_license_disclaimer(${file_name} ${CPM_PACKAGES})
  file(APPEND ${file_name}
    "\n-----\nThe following software may be included in this product: lua (License: https://www.lua.org/license.html)\n"
  )
  file(READ "${SCAFFOLDER_INCLUDE_DIR}/dualmc/LICENSE" license_TEXT)
  file(APPEND ${file_name}
    "\n-----\nThe following software may be included in this product: dualmc. This software contains the following license and notice below:\n\n${license_TEXT}\n"
  )
endfunction()

if(SCAFFOLD_BUILD_BINARY)
  message(STATUS "[Build Binary]")
  file(GLOB MAIN_SOURCES
    "${SCAFFOLDER_INCLUDE_DIR}/*/*.c"
    "${SCAFFOLDER_INCLUDE_DIR}/*/*.cpp"
    "${SCAFFOLDER_SOURCE_DIR}/utils.cpp"
    "${SCAFFOLDER_SOURCE_DIR}/Scaffolder.cpp"
    "${OpenCTM_SOURCE_DIR}/lib/*.c"
    "${OpenCTM_SOURCE_DIR}/lib/liblzma/*.c"
    "${vcglib_SOURCE_DIR}/wrap/ply/plylib.cpp"
  )

  add_executable(${PROJECT_NAME} ${MAIN_SOURCES} "${SCAFFOLDER_SOURCE_DIR}/Main.cpp")
  scaffolder_configure_target(${PROJECT_NAME})

  add_executable("${PROJECT_NAME}.SliceTest"
    ${MAIN_SOURCES}
    "${SCAFFOLDER_SOURCE_DIR}/SliceTest.cpp"
    "${SCAFFOLDER_SOURCE_DIR}/utils.cpp"
  )
  scaffolder_configure_target("${PROJECT_NAME}.SliceTest")
  include(GNUInstallDirs)
  install(TARGETS "${PROJECT_NAME}" DESTINATION bin COMPONENT scaffolder)
  install(TARGETS "${PROJECT_NAME}.SliceTest" DESTINATION bin COMPONENT scaffolder)
endif()

if(SCAFFOLD_BUILD_PYTHON)
  message(STATUS "[Build PyScaffolder]")
  find_package(Python COMPONENTS Interpreter Development)
  find_package(PyBind11 REQUIRED)

  file(GLOB PYTHON_SOURCES
    "${SCAFFOLDER_SOURCE_DIR}/python/*.cpp"
    "${SCAFFOLDER_SOURCE_DIR}/utils.cpp"
    "${SCAFFOLDER_SOURCE_DIR}/Scaffolder.cpp"
    "${OpenCTM_SOURCE_DIR}/lib/*.c"
    "${OpenCTM_SOURCE_DIR}/lib/liblzma/*.c"
    "${vcglib_SOURCE_DIR}/wrap/ply/plylib.cpp"
  )

  pybind11_add_module("Py${PROJECT_NAME}" ${PYTHON_SOURCES})
  target_include_directories("Py${PROJECT_NAME}" PRIVATE
    "${SCAFFOLDER_SOURCE_DIR}/python/"
    "${PyBind11_INCLUDE_DIR}"
  )
  scaffolder_configure_target("Py${PROJECT_NAME}")
  include(GNUInstallDirs)
  install(TARGETS "Py${PROJECT_NAME}"
    COMPONENT python
    DESTINATION "${Python_SITELIB}"
  )
endif()

write_third_party_license(${CMAKE_CURRENT_BINARY_DIR}/third_party.txt)