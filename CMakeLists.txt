cmake_minimum_required(VERSION 3.2)
project(hds-util VERSION 1.0.0)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose build type (Debug|Release|RelWithDebInfo|MinSizeRel)" FORCE)
endif()

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_FIND_LIBRARY_SUFFIXES ".a;${CMAKE_FIND_LIBRARY_SUFFIXES}") # Prefer libz.a when both are available

find_package(savvy REQUIRED)

add_executable(hds-util src/main.cpp src/merge_writer.hpp)
target_compile_definitions(hds-util PUBLIC -DVERSION="${PROJECT_VERSION}")
target_link_libraries(hds-util savvy)

install(TARGETS hds-util RUNTIME DESTINATION bin)
