cmake_minimum_required (VERSION 3.8)

project ("Scaffolder")

set(CGAL_DO_NOT_WARN_ABOUT_CMAKE_BUILD_TYPE ON)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# libigl
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         OFF)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           OFF)

message(STATUS "PROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}")

find_package(LIBIGL REQUIRED)
find_package(TBB REQUIRED)
find_package(VCG REQUIRED)

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
	"${PROJECT_SOURCE_DIR}/include/vcglib/*.h"
	"${PROJECT_SOURCE_DIR}/include/vcglib/*/*.h"
	"${PROJECT_SOURCE_DIR}/include/vcglib/*/*/*.h"
	"${PROJECT_SOURCE_DIR}/include/vcglib/*/*/*/*.h"
	"${VCG_INCLUDE_DIR}/wrap/ply/plylib.cpp"
)
include(update_deps_file)
update_deps_file("main_sources" "${MAIN_SOURCES}")

if(MSVC)
  # Enable parallel compilation for Visual Studio
  add_compile_options(/MP /bigobj)
  add_definitions(-DNOMINMAX -D_USE_MATH_DEFINES)
endif()

#add_definitions(-DDIP__IS_STATIC -DDIP__ENABLE_ASSERT -DDIP__HAS_JPEG -DDIP__EXCEPTIONS_RECORD_STACK_TRACE)
add_executable(${PROJECT_NAME} ${MAIN_SOURCES} ${PROJECT_SOURCE_DIR}/src/Main.cpp)
target_include_directories(${PROJECT_NAME} PRIVATE "${PROJECT_SOURCE_DIR}/include" ${TBB_INCLUDE_DIR} ${VCG_INCLUDE_DIR} ${DIPlib_INCLUDE_DIR} ${EIGEN_INCLUDE_DIR})
target_link_libraries(${PROJECT_NAME} PRIVATE igl::core DIP tbb_static)

add_executable(SliceTest ${MAIN_SOURCES} ${PROJECT_SOURCE_DIR}/src/SliceTest.cpp)
target_include_directories(SliceTest PRIVATE "${PROJECT_SOURCE_DIR}/include" ${TBB_INCLUDE_DIR} ${VCG_INCLUDE_DIR} ${EIGEN_INCLUDE_DIR})
target_link_libraries(SliceTest PRIVATE tbb_static)

add_executable(Fixer ${MAIN_SOURCES} ${PROJECT_SOURCE_DIR}/src/FixSelfIntersect.cpp)
target_include_directories(Fixer PRIVATE "${PROJECT_SOURCE_DIR}/include" ${TBB_INCLUDE_DIR} ${VCG_INCLUDE_DIR} ${EIGEN_INCLUDE_DIR})
target_link_libraries(Fixer tbb_static)
