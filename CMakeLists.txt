cmake_minimum_required(VERSION 3.13)
project(clustering)

set(CMAKE_CXX_STANDARD 20)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
#set(CMAKE_CXX_FLAGS_DEBUG "-O3")

set(ENV{OMP_NUM_THREADS} -1)


set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")
set (CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_LINKER_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")


#openMP stuff
#find_package(OpenMP REQUIRED)

#add_definitions(-D_GLIBCXX_PARALLEL)



find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
	link_libraries(OpenMP::OpenMP_CXX)
endif()

#find_package(OpenCV REQUIRED)
#include_directories( ${OpenCV_INCLUDE_DIRS} )

#set(CMAKE_CXX_FLAGS "-D\"HASVISUAL\"")

add_subdirectory(pybind11)
pybind11_add_module(klcluster src/geometry_types.cpp src/Curve.cpp src/io.cpp src/SparseFreespace.cpp src/Candidate.cpp src/center_clustering_algs.cpp src/python_wrapper.cpp)

add_executable(overhaul src/main.cpp src/geometry_types.cpp src/Curve.cpp src/CurveSimplification.cpp src/io.cpp src/SparseFreespace.cpp src/Candidate.cpp src/center_clustering_algs.cpp)
target_link_libraries(overhaul PUBLIC pybind11::embed)
#

#target_link_libraries(klcluster ${OpenCV_LIBS})
		#src/geometry_types.cpp  src/CurveSimplification.cpp src/io.cpp src/SparseFreespace.cpp src/FreespaceVisualizer.cpp src/Candidate.cpp src/center_clustering_algs.cpp src/python_wrapper.cpp)

#if(NOT TARGET OpenMP::OpenMP_CXX)
#	find_package(Threads REQUIRED)
#	add_library(OpenMP::OpenMP_CXX IMPORTED INTERFACE)
#	set_property(TARGET OpenMP::OpenMP_CXX
#			PROPERTY INTERFACE_COMPILE_OPTIONS ${OpenMP_CXX_FLAGS})
#	set_property(TARGET OpenMP::OpenMP_CXX
#			PROPERTY INTERFACE_LINK_LIBRARIES ${OpenMP_CXX_FLAGS} Threads::Threads)
#	link_libraries(OpenMP::OpenMP_CXX)
#endif()

#pybind11_add_module(backend src/main.cpp src/geometry_types.cpp src/Curve.cpp src/CurveSimplification.cpp src/io.cpp src/SparseFreespace.cpp src/FreespaceVisualizer.cpp src/Candidate.cpp src/center_clustering_algs.cpp)


#add_executable(overhaul src/main.cpp src/geometry_types.h src/basic_types.h src/id.h src/geometry_types.cpp src/Curve.cpp src/Curve.h src/CurveSimplification.cpp src/CurveSimplification.h src/io.cpp src/io.h src/defs.h src/SparseFreespace.cpp src/SparseFreespace.h src/FreespaceVisualizer.cpp src/FreespaceVisualizer.h src/Candidate.cpp src/Candidate.h src/center_clustering_algs.cpp src/center_clustering_algs.h src/SWatch.h)

#if(OpenCV_FOUND)
#	set(CMAKE_CXX_FLAGS "-D\"HASVISUAL\"")
#	include_directories( ${OpenCV_INCLUDE_DIRS} )
#	add_executable(overhaul src/main.cpp src/geometry_types.cpp src/Curve.cpp src/CurveSimplification.cpp src/io.cpp src/SparseFreespace.cpp src/FreespaceVisualizer.cpp src/Candidate.cpp src/center_clustering_algs.cpp)
#	target_link_libraries(overhaul ${OpenCV_LIBS})
#else()
#	set(CMAKE_CXX_FLAGS "-Wall")
#	set(CMAKE_CXX_FLAGS "-lstdc++")
#	add_executable(overhaul src/main.cpp src/geometry_types.h src/basic_types.h src/id.h src/geometry_types.cpp src/Curve.cpp src/Curve.h src/CurveSimplification.cpp src/CurveSimplification.h src/io.cpp src/io.h src/defs.h src/SparseFreespace.cpp src/SparseFreespace.h src/FreespaceVisualizer.cpp src/FreespaceVisualizer.h src/Candidate.cpp src/Candidate.h src/center_clustering_algs.cpp src/center_clustering_algs.h src/SWatch.h)
#endif()
