cmake_minimum_required(VERSION 3.2)
project(hds-util)

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_FIND_LIBRARY_SUFFIXES ".a;${CMAKE_FIND_LIBRARY_SUFFIXES}") # Prefer libz.a when both are available

find_package(savvy REQUIRED)

add_executable(hds-util src/main.cpp)
target_link_libraries(hds-util savvy)

install(TARGETS hds-util RUNTIME DESTINATION bin)
