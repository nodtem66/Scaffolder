cmake_minimum_required (VERSION 3.8)

project ("Scaffolder")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# libigl
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         OFF)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           OFF)
option(LIBIGL_USE_STATIC_LIBRARY "Use LibIGL as static library" ON)
option(IGL_STATIC_RUNTIME	"Use LibIGL as static library"	ON)
option(BUILD_PYSCAFFOLDER   "Build python library" ON)

if(MSVC)
  # Enable parallel compilation for Visual Studio
  add_compile_options(/MP /bigobj)
  add_definitions(-DNOMINMAX -D_USE_MATH_DEFINES -DWIN32_LEAN_AND_MEAN)
endif()

message(STATUS "CMAKE SYSTEM NAME: ${CMAKE_SYSTEM_NAME}")

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  message(STATUS "[OSX Detected]")	
  add_definitions(-D_DARWIN_C_SOURCE -D_GNU_SOURCE)
  add_compile_options(-D_XOPEN_SOURCE=600 -faligned-allocation)
endif()

if(${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
  message(STATUS "[CLANG Detected]")
endif()

message(STATUS "PROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}")
if(NOT EXISTS "${PROJECT_SOURCE_DIR}/external/libigl/CMakeLists.txt" OR
   NOT EXISTS "${PROJECT_SOURCE_DIR}/external/tbb/CMakeLists.txt" OR
   NOT EXISTS "${PROJECT_SOURCE_DIR}/external/sol2/CMakeLists.txt" OR
   NOT EXISTS "${PROJECT_SOURCE_DIR}/external/pybind11/CMakeLists.txt" OR
   NOT EXISTS "${PROJECT_SOURCE_DIR}/external/vcglib/wrap/")
  include(download_external_files)
  download_external_files(
    https://github.com/nodtem66/Scaffolder/releases/download/v1.5.0-alpha/external.tar.gz
	"${PROJECT_SOURCE_DIR}/scripts/external.tar.gz"
	e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855
  )
  file(ARCHIVE_EXTRACT INPUT "${PROJECT_SOURCE_DIR}/scripts/external.tar.gz" DESTINATION ${PROJECT_SOURCE_DIR})
endif()

find_package(Threads REQUIRED)
find_package(LIBIGL REQUIRED)
find_package(TBB REQUIRED)
find_package(VCG REQUIRED)
find_package(Sol2 REQUIRED)

# Add your project files
file(GLOB MAIN_SOURCES
	"${PROJECT_SOURCE_DIR}/include/*.c"
	"${PROJECT_SOURCE_DIR}/include/*.cpp"
	"${PROJECT_SOURCE_DIR}/include/*.tpp"
	"${PROJECT_SOURCE_DIR}/include/*.h"
	"${PROJECT_SOURCE_DIR}/include/*.hpp"
	"${PROJECT_SOURCE_DIR}/include/*/*.c"
	"${PROJECT_SOURCE_DIR}/include/*/*.cpp"
	"${PROJECT_SOURCE_DIR}/include/*/*.tpp"
	"${PROJECT_SOURCE_DIR}/include/*/*.h"
	"${PROJECT_SOURCE_DIR}/include/*/*.hpp"
	"${PROJECT_SOURCE_DIR}/src/utils.cpp"
	"${VCG_INCLUDE_DIR}/wrap/ply/plylib.cpp"
)
include(update_deps_file)
update_deps_file("main_sources" "${MAIN_SOURCES}")

# Build binary
add_executable(${PROJECT_NAME} ${MAIN_SOURCES} ${PROJECT_SOURCE_DIR}/src/Main.cpp)
target_include_directories(${PROJECT_NAME} PRIVATE "${PROJECT_SOURCE_DIR}/include" ${TBB_INCLUDE_DIR} ${VCG_INCLUDE_DIR} ${EIGEN_INCLUDE_DIR} ${SOL2_INCLUDE_DIR} ${LUA_INCLUDE_DIR})
target_link_libraries(${PROJECT_NAME} PRIVATE igl::core tbb Threads::Threads ${LUA_LIBRARIES})

add_executable("${PROJECT_NAME}.SliceTest" ${MAIN_SOURCES} ${PROJECT_SOURCE_DIR}/src/SliceTest.cpp)
target_include_directories("${PROJECT_NAME}.SliceTest" PRIVATE "${PROJECT_SOURCE_DIR}/include" ${VCG_INCLUDE_DIR} ${EIGEN_INCLUDE_DIR})
target_link_libraries("${PROJECT_NAME}.SliceTest" PRIVATE tbb)

# Install binary
include(GNUInstallDirs)
install(TARGETS "${PROJECT_NAME}" DESTINATION bin COMPONENT scaffolder)
install(TARGETS "${PROJECT_NAME}.SliceTest" DESTINATION bin COMPONENT scaffolder)

if(BUILD_PYSCAFFOLDER)
	message(STATUS "[Build PyScaffolder]")
	find_package(PyBind11 REQUIRED)

	file(GLOB PYTHON_SOURCES
		"${PROJECT_SOURCE_DIR}/src/python/*.cpp"
		"${PROJECT_SOURCE_DIR}/src/python/*.hpp"
		"${PROJECT_SOURCE_DIR}/src/utils.cpp"
	)
	include(update_deps_file)
	update_deps_file("python_sources" "${PYTHON_SOURCES}")

	pybind11_add_module("Py${PROJECT_NAME}" ${PYTHON_SOURCES})
	target_include_directories("Py${PROJECT_NAME}" PRIVATE "${PROJECT_SOURCE_DIR}/include" ${TBB_INCLUDE_DIR} ${VCG_INCLUDE_DIR} ${EIGEN_INCLUDE_DIR} ${SOL2_INCLUDE_DIR} ${LUA_INCLUDE_DIR} ${PyBind11_INCLUDE_DIR})
	target_link_libraries("Py${PROJECT_NAME}" PRIVATE igl::core tbb Threads::Threads ${LUA_LIBRARIES})
	install(TARGETS "Py${PROJECT_NAME}" DESTINATION lib)
endif()