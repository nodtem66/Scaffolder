cmake_minimum_required (VERSION 3.8)

project ("Scaffolder")

set(CGAL_DO_NOT_WARN_ABOUT_CMAKE_BUILD_TYPE ON)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# libigl
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         OFF)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           OFF)

message(STATUS "PROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}")

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
  add_definitions(-DNOMINMAX -D_USE_MATH_DEFINES -DWIN32_LEAN_AND_MEAN)
endif()

add_executable(${PROJECT_NAME} ${MAIN_SOURCES} ${PROJECT_SOURCE_DIR}/src/Main.cpp)
target_include_directories(${PROJECT_NAME} PRIVATE "${PROJECT_SOURCE_DIR}/include" ${TBB_INCLUDE_DIR} ${VCG_INCLUDE_DIR} ${EIGEN_INCLUDE_DIR} ${SOL2_INCLUDE_DIR} ${LUA_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} PRIVATE igl::core tbb Threads::Threads ${LUA_LIBRARIES})

#add_executable(Sol2_Test ${MAIN_SOURCES} ${PROJECT_SOURCE_DIR}/src/Sol2.cpp)
#target_include_directories(Sol2_Test PRIVATE "${PROJECT_SOURCE_DIR}/include" ${SOL2_INCLUDE_DIR} ${LUA_INCLUDE_DIRS})
#target_link_libraries(Sol2_Test PRIVATE Threads::Threads ${LUA_LIBRARIES})

add_executable(SliceTest ${MAIN_SOURCES} ${PROJECT_SOURCE_DIR}/src/SliceTest.cpp)
target_include_directories(SliceTest PRIVATE "${PROJECT_SOURCE_DIR}/include" ${VCG_INCLUDE_DIR} ${EIGEN_INCLUDE_DIR})
target_link_libraries(SliceTest PRIVATE tbb)

#add_executable(Fixer ${MAIN_SOURCES} ${PROJECT_SOURCE_DIR}/src/FixSelfIntersect.cpp)
#target_include_directories(Fixer PRIVATE "${PROJECT_SOURCE_DIR}/include" ${VCG_INCLUDE_DIR} ${EIGEN_INCLUDE_DIR})
#target_link_libraries(Fixer tbb)
