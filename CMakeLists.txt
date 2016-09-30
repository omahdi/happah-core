project(happah-core)
cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

# Libraries
find_package(GLM REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(Boost COMPONENTS iostreams REQUIRED)

# Paths & Files
file(GLOB_RECURSE SRC_FILES ${CMAKE_SOURCE_DIR}/lib/*.cpp)
file(GLOB_RECURSE INCLUDE_FILES ${CMAKE_SOURCE_DIR}/lib/*.h)
add_executable(${PROJECT_NAME} ${SRC_FILES} ${INCLUDE_FILES})

# Link & include libraries
include_directories(${CMAKE_SOURCE_DIR}/lib ${EIGEN3_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS} ${GLM_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} ${Boost_LIBRARIES})

target_compile_features(${PROJECT_NAME} PRIVATE cxx_generic_lambdas)
