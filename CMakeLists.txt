cmake_minimum_required(VERSION 3.13)
project(clustering)

set(CMAKE_CXX_STANDARD 20)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")
set (CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_LINKER_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")

set(ENV{OMP_NUM_THREADS} 12)

#openMP stuff
find_package(OpenMP REQUIRED)

add_definitions(-D_GLIBCXX_PARALLEL)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
	link_libraries(OpenMP::OpenMP_CXX)
endif()
if(NOT TARGET OpenMP::OpenMP_CXX)
	find_package(Threads REQUIRED)
	add_library(OpenMP::OpenMP_CXX IMPORTED INTERFACE)
	set_property(TARGET OpenMP::OpenMP_CXX
			PROPERTY INTERFACE_COMPILE_OPTIONS ${OpenMP_CXX_FLAGS})
	set_property(TARGET OpenMP::OpenMP_CXX
			PROPERTY INTERFACE_LINK_LIBRARIES ${OpenMP_CXX_FLAGS} Threads::Threads)

endif()
link_libraries(OpenMP::OpenMP_CXX)

set(OpenCV_DIR /Users/styx/GitPulls/opencv/install/lib/cmake/opencv4)
find_package(OpenCV)
if(OpenCV_FOUND)

	set(CMAKE_CXX_FLAGS "-D\"HASVISUAL\"")

	include_directories( ${OpenCV_INCLUDE_DIRS} )

	add_executable(clustering
			src/main.cpp
			#src/center_algs.cpp
			src/center_clustering_algs.cpp
			src/clustering_algs.cpp
			src/curve.cpp
			src/curve_simplification.cpp
			src/geometric_algs.cpp
			src/filter.cpp
			src/frechet_light.cpp
			src/frechet_matching.cpp
			src/geometry_basics.cpp
			src/io.cpp
			src/random.cpp
			src/freespace.cpp
			src/free_space_visualizer.cpp
			src/candidate.cpp
			src/SWatch.h)

	#target_link_options(clustering PRIVATE "/opt/homebrew/Cellar/gmp/6.2.1_1/lib/libgmp.a")
	#target_link_options(clustering PRIVATE "/opt/homebrew/Cellar/gmp/6.2.1_1/lib/libgmpxx.a")
	target_link_libraries(clustering ${OpenCV_LIBS} )
else()
	set(CMAKE_CXX_FLAGS "-Wall")
	set(CMAKE_CXX_FLAGS "-lstdc++")

	add_executable(clustering
			src/main.cpp
			#src/center_algs.cpp
			src/center_clustering_algs.cpp
			src/clustering_algs.cpp
			src/curve.cpp
			src/curve_simplification.cpp
			src/geometric_algs.cpp
			src/filter.cpp
			src/frechet_light.cpp
			src/frechet_matching.cpp
			src/geometry_basics.cpp
			src/io.cpp
			src/random.cpp
			src/freespace.cpp
			src/candidate.cpp
			src/SWatch.h)
endif()
