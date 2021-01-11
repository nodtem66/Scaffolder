# - Try to find the SOL2 library
# Once done this will define
#
#  SOL2_FOUND - system has SOL2
#  SOL2_INCLUDE_DIR - **the** SOL2 include directory
if(SOL2_FOUND)
    return()
endif()

find_path(SOL2_INCLUDE_DIR sol/sol.hpp
    HINTS
        ${SOL2_DIR}
        ENV SOL2_DIR
    PATHS
        ${PROJECT_SOURCE_DIR}/external/sol2
        ${PROJECT_SOURCE_DIR}/../external/sol2
        ${PROJECT_SOURCE_DIR}/../../external/sol2
    PATH_SUFFIXES include
    NO_DEFAULT_PATH
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Sol2
    "\nThePhD/sol2 not found --- You can download it using:\n\tgit clone https://github.com/ThePhD/sol2 ${CMAKE_SOURCE_DIR}/../sol2"
    SOL2_INCLUDE_DIR)
mark_as_advanced(SOL2_INCLUDE_DIR)

set(SOL2_LUA_VERSION "5.3.5" CACHE STRING "The version of Lua needed. Can be 5.1, 5.2, 5.3, 5.4, LuaJIT, or a more specific 3-part version number for a specifc Lua (e.g., 5.3.4 or luajit-2.0.5)")
set(SOL2_BUILD_LUA TRUE CACHE BOOL "Always build Lua, do not search for it in the system")
set(SOL2_PLATFORM "x64" CACHE STRING "Target platform to compile for when building binaries (x86, x64)")

# # # Platform
# Detect x86 and x64 stuff
if (SOL2_PLATFORM MATCHES "i686" OR SOL2_PLATFORM STREQUAL "x86")
	set(SOL2_IS_X86 TRUE)
elseif (SOL2_PLATFORM MATCHES "ARM64")
	set(IS_ARM64 TRUE)
	set(IS_X64 TRUE)
elseif (SOL2_PLATFORM MATCHES "ARM")
	set(IS_ARM TRUE)
elseif (SOL2_PLATFORM MATCHES "x86_64" OR SOL2_PLATFORM STREQUAL "x64")
	set(IS_X64 TRUE)
else()
	set(IS_X64 TRUE)
endif()

if (PROJECT_SOURCE_DIR STREQUAL ${CMAKE_SOURCE_DIR})
	set(SOL2_IS_TOP_LEVEL TRUE)
	message(STATUS "sol2 is the top-level directory...")
endif()

# # # sol2 Source Groups
# # Sources everyone is going to need
# Header files
file(GLOB SOL2_HEADER_SOURCES ${SOL2_INCLUDE_DIR}/sol*.hpp)
source_group(sol2 FILES ${SOL2_HEADER_SOURCES})

list(APPEND CMAKE_MODULE_PATH "${SOL2_INCLUDE_DIR}/../cmake/Modules")
message(STATUS "USE DIR: ${SOL2_INCLUDE_DIR}")

# # # Libraries
# Here, we pull in all the necessary libraries for building examples and tests
# Find threading library
if (NOT MSVC)
	if (SOL2_IS_X86)
		set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -m32")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m32")
		set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -m32")
		set(CMAKE_EXECUTABLE_LINKER_FLAGS "${CMAKE_EXECUTABLE_LINKER_FLAGS} -m32")
	endif()
	set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
	set(THREADS_PREFER_PTHREAD_FLAG TRUE)
else()
	string(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_C_FLAGS ${CMAKE_C_FLAGS})
	string(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS})
	if (BUILD_LUA_AS_DLL)
		string(REGEX REPLACE "/MT" "/MD" CMAKE_C_FLAGS ${CMAKE_C_FLAGS})
		string(REGEX REPLACE "/MT" "/MD" CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS})
	else ()
		string(REGEX REPLACE "/MD" "/MT" CMAKE_C_FLAGS ${CMAKE_C_FLAGS})
		string(REGEX REPLACE "/MD" "/MT" CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS})
	endif()
endif()
find_package(Threads REQUIRED)

string(TOLOWER ${SOL2_LUA_VERSION} NORMALIZED_LUA_VERSION)
# Find way to get Lua: build if requested, or attempt to build if no matching version is found
if (SOL2_BUILD_LUA)
	find_package(LuaBuild REQUIRED COMPONENTS ${SOL2_LUA_VERSION})
elseif (NOT SOL2_LUA_VERSION)
	find_package(LuaBuild REQUIRED)
else ()
	if (NORMALIZED_LUA_VERSION MATCHES "5.1")
		set(CREATE_LUALIB_TARGET TRUE)
		find_package(Lua 5.1 EXACT REQUIRED)
	elseif(NORMALIZED_LUA_VERSION MATCHES "5.2")
		set(CREATE_LUALIB_TARGET TRUE)
		find_package(Lua 5.2 EXACT REQUIRED)
	elseif(NORMALIZED_LUA_VERSION MATCHES "5.3")
		set(CREATE_LUALIB_TARGET TRUE)
		find_package(Lua 5.3 EXACT REQUIRED)
	elseif(NORMALIZED_LUA_VERSION MATCHES "5.4")
		set(CREATE_LUALIB_TARGET TRUE)
		find_package(Lua 5.4 EXACT REQUIRED)
	elseif(NORMALIZED_LUA_VERSION MATCHES "luajit")
		set(CREATE_LUALIB_TARGET TRUE)
		find_package(LuaJIT REQUIRED)
	else()
		find_package(LuaBuild ${SOL2_LUA_VERSION} REQUIRED)
	endif()
endif()

if (CREATE_LUALIB_TARGET AND LUA_FOUND)
	set(lualib lua_imported_lib_${SOL2_LUA_VERSION})
	foreach(lua_search_lib ${LUA_LIBRARIES})
		get_filename_component(lsl_fname ${lua_search_lib} NAME)
		if (lsl_fname MATCHES "lua" OR lsl_fname MATCHES "Lua" OR lsl_fname MATCHES "LUA")
			if (lsl_fname MATCHES "\.so|\.dylib|\.dll")
				set(lualibtype SHARED)
				set(lualiblocation ${lua_search_lib})
			else()
				set(lualibtype STATIC)
				set(lualiblocation ${lua_search_lib})
			endif()
		else()
			set(LUA_SEARCH_DEPENDENCY_LIBS ${LUA_SEARCH_DEPENDENCY_LIBS} "${lua_search_lib}")
		endif()
	endforeach()
	add_library(${lualib} ${lualibtype} IMPORTED)
	set_target_properties(${lualib}
		PROPERTIES 
		INTERFACE_INCLUDE_DIRECTORIES ${LUA_INCLUDE_DIR}
		INTERFACE_LINK_LIBRARIES ${LUA_SEARCH_DEPENDENCY_LIBS}
		IMPORTED_LINK_INTERFACE_LANGUAGES C
		IMPORTED_LOCATION ${lualiblocation})
	if (CMAKE_DL_LIBS)
		set_property(TARGET ${lualib}
			APPEND PROPERTY INTERFACE_LINK_LIBRARIES ${CMAKE_DL_LIBS})
	endif()
	set(LUA_LIBRARIES ${lualib})
endif()

if (NOT LUA_FOUND AND NOT LUABUILD_FOUND)
	message(FATAL_ERROR "sol2 Lua \"${SOL2_LUA_VERSION}\" not found and could not be targeted for building")
endif()
message(STATUS "USE DIR: ${LUA_LIBRARIES} ${LUA_INCLUDE_DIRS}")
