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

message(STATUS "PROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}")
if(NOT EXISTS "${PROJECT_SOURCE_DIR}/external/libigl/CMakeLists.txt" OR
   NOT EXISTS "${PROJECT_SOURCE_DIR}/external/tbb/CMakeLists.txt" OR
   NOT EXISTS "${PROJECT_SOURCE_DIR}/external/sol2/CMakeLists.txt" OR
   NOT EXISTS "${PROJECT_SOURCE_DIR}/external/vcglib/wrap/")
  include(download_external_files)
  download_external_files(
    https://github.com/nodtem66/Scaffolder/releases/download/v1.5.0-alpha/external.tar.gz
	"${PROJECT_SOURCE_DIR}/scripts/external.tar.gz"
	40575b7a8486ce66e2f2551cb94fff11d9f590fb5924d50db4b455a3bc1e299f
  )
  file(ARCHIVE_EXTRACT INPUT "${PROJECT_SOURCE_DIR}/scripts/external.tar.gz" DESTINATION ${PROJECT_SOURCE_DIR})
endif()

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
	"${VCG_INCLUDE_DIR}/wrap/ply/plylib.cpp"
)
include(update_deps_file)
update_deps_file("main_sources" "${MAIN_SOURCES}")

if(MSVC)
  # Enable parallel compilation for Visual Studio
  add_compile_options(/MP /bigobj)
  add_definitions(-DNOMINMAX -D_USE_MATH_DEFINES -DWIN32_LEAN_AND_MEAN)
endif()

add_executable(${PROJECT_NAME} ${MAIN_SOURCES} ${PROJECT_SOURCE_DIR}/src/Main.cpp )
target_include_directories(${PROJECT_NAME} PRIVATE "${PROJECT_SOURCE_DIR}/include" ${TBB_INCLUDE_DIR} ${VCG_INCLUDE_DIR} ${EIGEN_INCLUDE_DIR} ${SOL2_INCLUDE_DIR} ${LUA_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} PRIVATE igl::core tbb Threads::Threads ${LUA_LIBRARIES})

#add_executable(Sol2_Test ${MAIN_SOURCES} ${PROJECT_SOURCE_DIR}/src/Sol2.cpp)
#target_include_directories(Sol2_Test PRIVATE "${PROJECT_SOURCE_DIR}/include" ${SOL2_INCLUDE_DIR} ${LUA_INCLUDE_DIRS})
#target_link_libraries(Sol2_Test PRIVATE Threads::Threads ${LUA_LIBRARIES})

add_executable("${PROJECT_NAME}.SliceTest" ${MAIN_SOURCES} ${PROJECT_SOURCE_DIR}/src/SliceTest.cpp )
target_include_directories("${PROJECT_NAME}.SliceTest" PRIVATE "${PROJECT_SOURCE_DIR}/include" ${VCG_INCLUDE_DIR} ${EIGEN_INCLUDE_DIR})
target_link_libraries("${PROJECT_NAME}.SliceTest" PRIVATE tbb)

#add_executable(Fixer ${MAIN_SOURCES} ${PROJECT_SOURCE_DIR}/src/FixSelfIntersect.cpp)
#target_include_directories(Fixer PRIVATE "${PROJECT_SOURCE_DIR}/include" ${VCG_INCLUDE_DIR} ${EIGEN_INCLUDE_DIR})
#target_link_libraries(Fixer tbb)

include(GNUInstallDirs)
install(TARGETS "${PROJECT_NAME}" DESTINATION bin COMPONENT scaffolder)
install(TARGETS "${PROJECT_NAME}.SliceTest" DESTINATION bin COMPONENT scaffolder)