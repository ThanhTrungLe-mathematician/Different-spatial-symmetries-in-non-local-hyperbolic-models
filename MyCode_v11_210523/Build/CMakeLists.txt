cmake_minimum_required(VERSION 2.8.12)
project (Main)                  #name of project

#parallel computing
OPTION (USE_OpenMP "Use OpenMP" ON)
IF(USE_OpenMP)
    FIND_PACKAGE(OpenMP)
IF(OPENMP_FOUND)
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
ENDIF()
ENDIF()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp")


include_directories(../Source/Header)     #folder contains all .h files
file(GLOB SOURCES "../Source/Src/*.cpp")  #folder contains all .cpp files responding .h files
#set(SOURCES ../Src/utilities.cpp)
set(MAIN ../Main/Main.cpp)         #set file main.cpp 
add_executable(Main ${MAIN} ${SOURCES})

add_custom_target(run
            COMMAND ./Main
            DEPENDS Main
            )

add_custom_target(clean_cmake
            COMMAND rm -f Main
            COMMAND rm -f CMakeCache.txt
            COMMAND rm -f cmake_install.cmake
            COMMAND rm -rf CMakeFiles
            COMMAND rm -f Makefile
            )

add_custom_target(clean_txt
            find ../Output/ -name '*.txt' -type f -delete
            )

add_custom_target(clean_nohup
            COMMAND rm -f nohup.out
            )
