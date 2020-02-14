cmake_minimum_required (VERSION 3.8)

project ("Scaffolder")

set(CGAL_DO_NOT_WARN_ABOUT_CMAKE_BUILD_TYPE ON)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# libigl
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)

find_package(LIBIGL REQUIRED)
find_package(VCG REQUIRED)
find_package(DIPlib REQUIRED)
#find_package(HDF5 CONFIG REQUIRED)
#find_package(HighFive REQUIRED)

# Add your project files
file(GLOB SRCFILES src/*.cpp)
add_definitions(-DNOMINMAX -D_USE_MATH_DEFINES -DDIP__IS_STATIC -DDIP__ENABLE_ASSERT -DDIP__HAS_JPEG -DDIP__EXCEPTIONS_RECORD_STACK_TRACE)
add_executable(${PROJECT_NAME} ${SRCFILES} ${VCG_INCLUDE_DIR}/wrap/ply/plylib.cpp)
target_include_directories(${PROJECT_NAME} PRIVATE "include" ${VCG_INCLUDE_DIR} ${DIPlib_INCLUDE_DIR})
target_link_libraries(${PROJECT_NAME} PRIVATE igl::core DIP)
