cmake_minimum_required (VERSION 3.8)

project ("Scaffolder")

set(CGAL_DO_NOT_WARN_ABOUT_CMAKE_BUILD_TYPE ON)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# libigl
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         OFF)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           OFF)

find_package(LIBIGL REQUIRED)
find_package(VCG REQUIRED)
find_package(DIPlib REQUIRED)
find_package(TBB REQUIRED)
#find_package(HDF5 CONFIG REQUIRED)
#find_package(HighFive REQUIRED)
#find_package(OpenMP)
#if (OpenMP_FOUND)
#    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
#	if (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
#		string(REPLACE "-openmp" "-openmp:experimental" CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS})
#	endif()
#endif()

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
)
include(update_deps_file)
update_deps_file("main_sources" "${MAIN_SOURCES}")

add_definitions(-DNOMINMAX -D_USE_MATH_DEFINES -DDIP__IS_STATIC -DDIP__ENABLE_ASSERT -DDIP__HAS_JPEG -DDIP__EXCEPTIONS_RECORD_STACK_TRACE)
add_executable(${PROJECT_NAME} ${MAIN_SOURCES} ${PROJECT_SOURCE_DIR}/src/Main.cpp ${VCG_INCLUDE_DIR}/wrap/ply/plylib.cpp)
target_include_directories(${PROJECT_NAME} PRIVATE "${PROJECT_SOURCE_DIR}/include" ${VCG_INCLUDE_DIR} ${DIPlib_INCLUDE_DIR})
target_link_libraries(${PROJECT_NAME} PRIVATE igl::core DIP)

add_executable(SliceTest ${MAIN_SOURCES} ${PROJECT_SOURCE_DIR}/src/SliceTest.cpp)
target_include_directories(SliceTest PRIVATE "${PROJECT_SOURCE_DIR}/include" ${VCG_INCLUDE_DIR} ${TBB_INCLUDE_DIR})
target_link_libraries(SliceTest PRIVATE igl::core tbb_static)
