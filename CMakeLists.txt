cmake_minimum_required(VERSION 3.12)

find_package(MPI REQUIRED)
find_package(Boost COMPONENTS
    mpi
    iostreams
    serialization
    program_options
    REQUIRED)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(mpi-pi-lib
    "include/common.hpp"
	"include/random.hpp"
	"include/pi.hpp"
	"src/common.cpp")

target_include_directories(mpi-pi-lib PUBLIC
    ${Boost_INCLUDE_DIRS}
    ${MPI_INCLUDE_PATH})
target_link_libraries(mpi-pi-lib PRIVATE
    ${Boost_LIBRARIES}
    ${MPI_LIBRARIES})

add_executable(mpi-pi
	"config.hpp"
	"main.cpp")

target_include_directories(mpi-pi PRIVATE
	${Boost_INCLUDE_DIRS}
	${MPI_INCLUDE_PATH})

target_link_libraries(mpi-pi PRIVATE
	mpi-pi-lib
	${Boost_LIBRARIES}
	${MPI_LIBRARIES})
