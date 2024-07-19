cmake_minimum_required(VERSION 3.15...3.27)

set(CMAKE_CXX_STANDARD 17 CACHE STRING "The C++ standard to use")
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

project(FLOOD VERSION 1.0 DESCRIPTION "synthetic car sequence recorder" LANGUAGES CXX)

# TARGET
set(FLOOD_SOURCE_FILES
    "src/main.cpp"
)
add_executable(FLOOD ${FLOOD_SOURCE_FILES})
target_include_directories(FLOOD PUBLIC "src")

set(RAYLIB_VERSION 5.0)
find_package(raylib ${RAYLIB_VERSION} QUIET) # QUIET or REQUIRED
if (NOT raylib_FOUND) # If there's none, fetch and build raylib
  include(FetchContent)
  FetchContent_Declare(
    raylib
    DOWNLOAD_EXTRACT_TIMESTAMP OFF
    URL https://github.com/raysan5/raylib/archive/refs/tags/${RAYLIB_VERSION}.tar.gz
  )
  FetchContent_GetProperties(raylib)
  if (NOT raylib_POPULATED) # Have we downloaded raylib yet?
    set(FETCHCONTENT_QUIET NO)
    FetchContent_Populate(raylib)
    set(BUILD_EXAMPLES OFF CACHE BOOL "" FORCE) # don't build the supplied examples
    add_subdirectory(${raylib_SOURCE_DIR} ${raylib_BINARY_DIR})
  endif()
endif()
target_link_libraries(FLOOD PUBLIC raylib)

add_custom_target(prepare_rel
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/Release
)
add_dependencies(raylib prepare_rel)

add_custom_target(copy_assets_rel
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_LIST_DIR}/res ${CMAKE_CURRENT_BINARY_DIR}/Release/res
)
add_dependencies(raylib copy_assets_rel)

add_custom_target(prepare_deb
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/Debug
)
add_dependencies(raylib prepare_deb)

add_custom_target(copy_assets_deb
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_LIST_DIR}/res ${CMAKE_CURRENT_BINARY_DIR}/Debug/res
)
add_dependencies(raylib copy_assets_deb)