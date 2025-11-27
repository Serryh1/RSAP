cmake_minimum_required(VERSION 3.5.0)
project(RSAP VERSION 0.1.0 LANGUAGES C CXX)

# if(${CMAKE_BUILD_TYPE} STREQUAL "Debug")
#     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
# endif()
set(CMAKE_CXX_STANDARD 17)

# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address,undefined")
# set(CMAKE_BUILD_TYPE Release)

add_executable(RSAP SRC/main.cpp SRC/GW.cpp  SRC/SF-Solution.cpp SRC/SF-Functions.cpp  SRC/Input.cpp SRC/BipartiteGraph.cpp SRC/HRA.cpp)


