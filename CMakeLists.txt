cmake_minimum_required(VERSION 3.14)
project(AMES)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Prefer the active environment/PATH over a newer system Python.
# An explicit -DPython_EXECUTABLE=/path/to/python takes precedence.
set(Python_FIND_STRATEGY LOCATION)
set(Python_FIND_FRAMEWORK LAST)
find_package(Python REQUIRED COMPONENTS Interpreter Development)
message(STATUS "AMES Python: ${Python_EXECUTABLE} (${Python_VERSION})")
message(STATUS "AMES Python install directory: ${Python_SITEARCH}")

# Find OpenMP (required for C/C++)

# OpenMP Configuration (macOS-specific)
if(APPLE)
    # Share the selected Python environment's runtime (e.g. NumPy/Conda).
    # Allow an explicit -DAMES_OPENMP_ROOT=/path/to/prefix override.
    execute_process(COMMAND "${Python_EXECUTABLE}" -c "import sys; print(sys.prefix)"
        OUTPUT_VARIABLE AMES_PYTHON_PREFIX OUTPUT_STRIP_TRAILING_WHITESPACE
        RESULT_VARIABLE AMES_PYTHON_PREFIX_STATUS)
    if(NOT AMES_PYTHON_PREFIX_STATUS EQUAL 0)
        message(FATAL_ERROR "Unable to query the selected Python environment")
    endif()
    if(NOT AMES_OPENMP_ROOT)
        if(EXISTS "${AMES_PYTHON_PREFIX}/lib/libomp.dylib" AND
           EXISTS "${AMES_PYTHON_PREFIX}/include/omp.h")
            set(AMES_OPENMP_ROOT "${AMES_PYTHON_PREFIX}")
        else()
            set(AMES_OPENMP_ROOT "/opt/homebrew/opt/libomp")
        endif()
    endif()
    message(STATUS "AMES OpenMP runtime: ${AMES_OPENMP_ROOT}/lib/libomp.dylib")
    set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp -I${AMES_OPENMP_ROOT}/include")
    set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I${AMES_OPENMP_ROOT}/include")
    set(OpenMP_C_LIB_NAMES "omp")
    set(OpenMP_CXX_LIB_NAMES "omp")
    set(OpenMP_omp_LIBRARY "${AMES_OPENMP_ROOT}/lib/libomp.dylib" CACHE FILEPATH "AMES OpenMP runtime" FORCE)
    
    find_package(OpenMP REQUIRED)
    
    # Special handling for AppleClang
    if(CMAKE_CXX_COMPILER_ID MATCHES "AppleClang")
        add_compile_options(-Xpreprocessor)
    endif()
else()
    # Standard OpenMP for other platforms
    find_package(OpenMP REQUIRED)
endif()


set(PYBIND11_FINDPYTHON ON)
find_package(pybind11 REQUIRED)

# Add your executable/library

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -ffast-math")

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DNDEBUG")

# 3. Configure OpenMP specifically for PyBind11
if(OpenMP_CXX_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    add_definitions(-DUSE_OPENMP)
endif()

set_property(SOURCE AMES.i PROPERTY CPLUSPLUS ON)

pybind11_add_module(${PROJECT_NAME}
src/module/main_bindings.cpp
src/module/utility.cpp
src/module/particle.cpp
src/module/source.cpp
src/module/synchrotron.cpp
src/module/compton.cpp
src/module/gammagamma.cpp
src/module/electrondistribution.cpp
src/module/photonbackground.cpp
src/GRBAfterglow.cpp
)

target_include_directories(AMES PRIVATE 
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

if(APPLE)
    # Check if using Homebrew LLVM Clang
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        target_compile_options(${PROJECT_NAME} PRIVATE
            $<$<COMPILE_LANGUAGE:CXX>:
                -Xpreprocessor -fopenmp
                -fopenmp-simd
                -ffast-math
                -fvisibility=default
                -fno-trapping-math
                -O3
                -march=native
                # Uncomment these for vectorization diagnostics (Clang-only)
                #-Rpass=loop-vectorize
                #-Rpass-missed=loop-vectorize
                #-Rpass-analysis=loop-vectorize
            >
        )
    endif()
else()  # Linux or other UNIX
    target_compile_options(${PROJECT_NAME} PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:
            -fopenmp
            -fopenmp-simd
            -ffast-math
            -fvisibility=default
            -fno-trapping-math
            -O3
            -march=native
        >
    )
endif()

# Conditional compilation for Photopair
target_compile_definitions(${PROJECT_NAME} PRIVATE
    -DPHOTOPAIR_USE_OPENMP
)

set_property(TARGET ${PROJECT_NAME} PROPERTY SWIG_INCLUDE_DIRECTORIES ${CMAKE_CURRENT_SOURCE_DIR}/include)
set_property(TARGET ${PROJECT_NAME} PROPERTY INCLUDE_DIRECTORIES ${CMAKE_CURRENT_SOURCE_DIR}/include ${Python_INCLUDE_DIRS})

add_definitions(-DDATA_PATH=\"${CMAKE_CURRENT_SOURCE_DIR}\" -DINFO)
#target_link_libraries(${PROJECT_NAME} PRIVATE Python::Python)
#target_link_libraries(${PROJECT_NAME} PRIVATE Python3::Python GSL::gsl GSL::gslcblas)
#target_link_libraries(${PROJECT_NAME} PRIVATE ${Python_LIBRARIES})
# pybind11_add_module supplies Python extension linking. Linking libpython
# explicitly can load a second Python runtime with Conda on macOS.
target_link_libraries(${PROJECT_NAME} PRIVATE OpenMP::OpenMP_CXX) #syu

#target_link_libraries(${PROJECT_NAME} PRIVATE ${Python_LIBRARIES} GSL::gsl GSL::gslcblas)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
  target_link_libraries(${PROJECT_NAME} ${Boost_LIBRARIES})
endif()

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/include  # For your headers
    ${Python_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIRS}
)

set(CMAKE_CXX_VISIBILITY_PRESET "default")
set(CMAKE_VISIBILITY_INLINES_HIDDEN OFF)
if(APPLE)
    # we change the suffix because swig generates a file expecting a .so on mac.
#    set_target_properties(${PROJECT_NAME} PROPERTIES SUFFIX ".so")
    set_target_properties(${PROJECT_NAME} PROPERTIES CXX_VISIBILITY_PRESET "default")
endif()

install(TARGETS ${PROJECT_NAME}
    LIBRARY DESTINATION ${Python_SITEARCH}
)
