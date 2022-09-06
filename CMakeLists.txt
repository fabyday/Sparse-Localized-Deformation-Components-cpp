cmake_minimum_required(VERSION 3.16)

project(Sparse-Localized-Deformation-Components-cpp)
cmake_policy(SET CMP0074 NEW)


#if(CMAKE_CONFIGURATION_TYPES)
#	if(NOT "DEBUG_WITH_LOG" IN_LIST CMAKE_CONFIGURATION_TYPES)
#		list(APPEND CMAKE_CONFIGURATION_TYPES "DEBUG_WITH_LOG")
#	endif()
#endif()



list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmakes)



add_compile_definitions("_CRT_SECURE_NO_WARNINGS")
add_compile_definitions("WIN32")
add_compile_definitions("_WINDOWS")
add_compile_definitions("IMGUI_IMPL_OPENGL_LOADER_GLAD")

set(CMAKE_CXX_STANDARD 17)


# Libigl
option(LIBIGL_WITH_COMISO            "Use CoMiso"                   ON)
option(LIBIGL_WITH_EMBREE            "Use Embree"                   ON)
option(LIBIGL_WITH_OPENGL            "Use OpenGL"                   ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"                     ON)
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI "Use ImGui"                    ON)
option(LIBIGL_WITH_PNG               "Use PNG"                      ON)
option(LIBIGL_WITH_TETGEN            "Use Tetgen"                   ON)
option(LIBIGL_WITH_TRIANGLE          "Use Triangle"                 ON)
option(LIBIGL_WITH_PREDICATES        "Use exact predicates"         ON)
option(LIBIGL_WITH_XML               "Use XML"                      ON)
option(PRINT_PRETTY_LOG_ALL               "Use matrix log"                      OFF)
option(PRINT_PRETTY_LOG_DENSE               "Use matrix log"                      OFF)



# library
include(libigl)
include(spdlog)
include(yaml-cpp)
#include(pcl)

#include(zlib)
#set(ZLIB_ROOT ${zlib_SOURCE_DIR})
#include(cnpy)

#use multi byte compile
add_definitions(-D_MBCS)

set(cpp_files src/splocs.cpp)

set(EXAMPLE_PATH "${CMAKE_SOURCE_DIR}/examples")

add_compile_definitions(
	EXAMPLE_PATH="${EXAMPLE_PATH}/"
	face205="${EXAMPLE_PATH}/face-poses/"
	navitan="${EXAMPLE_PATH}/lowpoly/cat/animation/"
	big_face="${EXAMPLE_PATH}/big-face/"
)

# set(CGAL_SOURCE_DIR "${CMAKE_SOURCE_DIR}/ThirdParty/CGAL-5.5")
# set(CGAL_DIR "${CMAKE_BINARY_DIR}/CGAL")
# add_subdirectory("ThirdParty/CGAL-5.5") 



#include ("${CMAKE_SOURCE_DIR}/ThirdParty/CGAL-5.5/lib/cmake/CGAL/CGALConfig.cmake")
#link_libraries(${CGAL_LIBRARIES})
#message(${CGAL_LIBRARIES})
#message(${CGAL_INCLUDE_DIRS})
#include_directories(${CGAL_INCLUDE_DIRS})




add_library(splocs ${cpp_files} "src/splocs.h" "src/splocs.cpp")
target_include_directories(splocs PUBLIC ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(splocs igl::core igl::opengl_glfw igl::opengl_glfw_imgui)


# excutable
add_subdirectory(tests)

