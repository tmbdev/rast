#CMake minimum version
cmake_minimum_required(VERSION 2.8)

#Project name
#--------------------------------------------
#Remember: Need to be changed in Project->settings->General->Executable to run
#          $(ProjectPath)build\sample
#--------------------------------------------
set(project_name RAST_test)


#Project definition
#--------------------------------------------
project( ${project_name} )

set(IULIB_INCLUDE_DIRS /home/shawnle/projects/iulib_git/)
set(IMGLIB_INCLUDE_DIRS /home/shawnle/projects/iulib_git/imglib)
set(IMGBITS_INCLUDE_DIRS /home/shawnle/projects/iulib_git/imgbits)
set(IMGIO_INCLUDE_DIRS /home/shawnle/projects/iulib_git/imgio)
set(COLIB_INCLUDE_DIRS /home/shawnle/projects/iulib_git/colib)
set(UTILS_INCLUDE_DIRS /home/shawnle/projects/iulib_git/utils)
set(VIDIO_INCLUDE_DIRS /home/shawnle/projects/iulib_git/vidio)
set(VIDIOUTILS_INCLUDE_DIRS /home/shawnle/projects/iulib_git/vidio/utils)
set(COMPONENTS_INCLUDE_DIRS /home/shawnle/projects/iulib_git/components)

#Find libraries
#--------------------------------------------
find_package( OpenCV REQUIRED )
#find_package( Boost REQUIRED )
#find_package( Eigen3 REQUIRED )
#find_package( OpenNI REQUIRED )
#find_package( Nite REQUIRED )


#Include directories
#--------------------------------------------
include_directories(${OpenCV_INCLUDE_DIRS})
#include_directories(${Boost_INCLUDE_DIRS})
#include_directories(${EIGEN3_INCLUDE_DIR})
#include_directories(${OpenNI_INCLUDE_DIR})
#include_directories(${Nite_INCLUDE_DIR})

## LOTS OF IULIB MODULES
include_directories(${COMPONENTS_INCLUDE_DIRS})
include_directories(${IULIB_INCLUDE_DIRS})
include_directories(${IMGLIB_INCLUDE_DIRS})
include_directories(${IMGBITS_INCLUDE_DIRS})
include_directories(${COLIB_INCLUDE_DIRS})
include_directories(${IMGIO_INCLUDE_DIRS})
include_directories(${UTILS_INCLUDE_DIRS})
include_directories(${VIDIO_INCLUDE_DIRS})
include_directories(${VIDIOUTILS_INCLUDE_DIRS})


#Link directories
#--------------------------------------------
link_directories(${Boost_LIBRARY_DIR})

#message(${OpenCV_INCLUDE_DIRS})
#message(${OpenCV_LIBS})


#Files to be compiled
#--------------------------------------------
#add_executable( ${project_name} calignmentp2d.cc rast-test.cc crastrs2d.cc crasts2d.cc cliness2d.cc clinesp2d.cc cinstancep2d.cc cedges.cc crastss2d.cc crastp2d.cc rast.cc) 
add_executable( ${project_name} rast-test.cc crasts2d.cc clinesp2d.cc cliness2d.cc crastp2d.cc cinstancep2d.cc EdgeDetector.cpp)


#Link libraries
#--------------------------------------------
target_link_libraries( ${project_name} ${OpenCV_LIBS})
#target_link_libraries( ${project_name} ${Boost_LIBRARIES})
#target_link_libraries( ${project_name} ${OpenNI_LIBRARIES})
#target_link_libraries( ${project_name} ${Nite_LIBRARIES})
#target_link_libraries( ${project_name} -fopenmp )
