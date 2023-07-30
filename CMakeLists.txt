cmake_minimum_required(VERSION 3.15)
project(AMES C CXX)

set(CMAKE_CXX_STANDARD 14)

find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})

find_package(Python COMPONENTS Interpreter Development)
#find_package(Python COMPONENTS Interpreter Development REQUIRED)
#find_package(Python 3.10 EXACT COMPONENTS Interpreter Development REQUIRED)

message(STATUS "Python_FOUND: ${Python_FOUND}")
message(STATUS "Python_EXECUTABLE: ${Python_EXECUTABLE}")
message(STATUS "Python_VERSION: ${Python_VERSION}")
message(STATUS "Python_INCLUDE_DIRS: ${Python_INCLUDE_DIRS}")
message(STATUS "Python_LIBRARY_DIRS: ${Python_LIBRARY_DIRS}")
message(STATUS "Python_LIBRARIES: ${Python_LIBRARIES}")

include_directories(${Python_INCLUDE_DIRS})

set_property(SOURCE AMES.i PROPERTY CPLUSPLUS ON) 

# Set the name of your SWIG module and the source files
set(SWIG_MODULE_myproject_NAME AMES)
set(SWIG_MODULE_myproject_EXTRA_DEPS ${CMAKE_CURRENT_SOURCE_DIR}/include)

swig_add_library(${SWIG_MODULE_myproject_NAME} 
    TYPE MODULE 
    LANGUAGE python 
    SOURCES AMES.i 
src/module/utility.cpp
src/module/particle.cpp
src/module/source.cpp
src/module/synchrotron.cpp
src/module/inversecompton.cpp
src/module/gammagamma.cpp
src/module/electrondistribution.cpp
src/module/photonbackground.cpp
src/GRBAfterglow.cpp
)



# Set include directories
set_property(TARGET ${SWIG_MODULE_myproject_NAME} PROPERTY SWIG_INCLUDE_DIRECTORIES ${CMAKE_CURRENT_SOURCE_DIR}/include)
set_property(TARGET ${SWIG_MODULE_myproject_NAME} PROPERTY INCLUDE_DIRECTORIES ${CMAKE_CURRENT_SOURCE_DIR}/include ${Python_INCLUDE_DIRS})

# Add any other definitions or flags you need
add_definitions(-DDATA_PATH=\"${CMAKE_CURRENT_SOURCE_DIR}\" -DINFO)

# Link the SWIG module with the necessary libraries
target_link_libraries(${SWIG_MODULE_myproject_NAME} PRIVATE ${Python_LIBRARIES})

# Install the SWIG module to the virtual environment
set(CMAKE_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
install(TARGETS ${SWIG_MODULE_myproject_NAME}
    DESTINATION ${CMAKE_INSTALL_PREFIX})

#if(APPLE)
    # we change the suffix because swig generates a file expecting a .so on mac.
#    set_target_properties(${SWIG_MODULE_myproject_NAME} PROPERTIES SUFFIX ".so")
#endif()

