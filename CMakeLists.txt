cmake_minimum_required(VERSION 3.14)
project(AMES)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#find_package(Python 3.13 COMPONENTS Interpreter Development)
find_package(Python COMPONENTS Interpreter Development)

# Find OpenMP (required for C/C++)

# OpenMP Configuration (macOS-specific)
if(APPLE)
    # Manually set OpenMP paths for Homebrew-installed libomp
    set(OPENMP_ROOT "/opt/homebrew/opt/libomp")
    set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp -I${OPENMP_ROOT}/include")
    set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I${OPENMP_ROOT}/include")
    set(OpenMP_C_LIB_NAMES "omp")
    set(OpenMP_CXX_LIB_NAMES "omp")
    set(OpenMP_omp_LIBRARY "${OPENMP_ROOT}/lib/libomp.dylib")
    
    find_package(OpenMP REQUIRED)
    
    # Special handling for AppleClang
    if(CMAKE_CXX_COMPILER_ID MATCHES "AppleClang")
        add_compile_options(-Xpreprocessor)
        link_directories(${OPENMP_ROOT}/lib)
    endif()
else()
    # Standard OpenMP for other platforms
    find_package(OpenMP REQUIRED)
endif()


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
        target_link_libraries(${PROJECT_NAME} PRIVATE omp)
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
target_link_libraries(${PROJECT_NAME} PRIVATE ${Python_LIBRARIES})
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
