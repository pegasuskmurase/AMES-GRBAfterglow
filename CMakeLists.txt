cmake_minimum_required(VERSION 3.18)
project(AMES LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Preserve an optimized default without forcing -O3 into Debug builds.
if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
        Debug Release RelWithDebInfo MinSizeRel)
endif()

option(AMES_NATIVE_ARCH
    "Optimize for the build CPU; do not use for portable or universal binaries" OFF)
option(AMES_FAST_MATH
    "Enable the historical fast-math setting; may affect floating-point results" ON)
set(AMES_OPENMP_ROOT "" CACHE PATH
    "Optional macOS Clang OpenMP prefix containing include/omp.h and lib/libomp.dylib")

# Keep the chosen interpreter consistent with the Python headers and install path.
# An explicit -DPython_EXECUTABLE=/absolute/path/to/python takes precedence.
set(PYBIND11_FINDPYTHON ON)
set(Python_FIND_STRATEGY LOCATION)
set(Python_FIND_FRAMEWORK LAST)
find_package(Python 3 REQUIRED COMPONENTS Interpreter Development.Module)

execute_process(
    COMMAND "${Python_EXECUTABLE}" -c "import sys; print(sys.prefix)"
    OUTPUT_VARIABLE AMES_PYTHON_PREFIX
    OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE _ames_python_prefix_status)
execute_process(
    COMMAND "${Python_EXECUTABLE}" -c "import platform; print(platform.machine())"
    OUTPUT_VARIABLE AMES_PYTHON_ARCH
    OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE _ames_python_arch_status)
if(NOT _ames_python_prefix_status EQUAL 0 OR
   NOT _ames_python_arch_status EQUAL 0)
    message(FATAL_ERROR "Unable to query the selected Python interpreter.")
endif()
if(NOT Python_SITEARCH)
    message(FATAL_ERROR "The selected Python did not report a site-packages directory.")
endif()

message(STATUS "AMES Python: ${Python_EXECUTABLE} (${Python_VERSION})")
message(STATUS "AMES Python process architecture: ${AMES_PYTHON_ARCH}")
message(STATUS "AMES Python install directory: ${Python_SITEARCH}")
message(STATUS "AMES C++ compiler: ${CMAKE_CXX_COMPILER} (${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION})")

if(APPLE)
    if(CMAKE_OSX_ARCHITECTURES)
        message(STATUS "AMES macOS target architectures: ${CMAKE_OSX_ARCHITECTURES}")
        if(NOT AMES_PYTHON_ARCH IN_LIST CMAKE_OSX_ARCHITECTURES)
            message(WARNING
                "Python is running as ${AMES_PYTHON_ARCH}, while the requested "
                "macOS target architectures are ${CMAKE_OSX_ARCHITECTURES}. "
                "For a native installation, select a matching Python and configure "
                "a fresh build directory. A universal interpreter can run as only "
                "one of its supported slices, so this check is diagnostic only.")
        endif()
    else()
        message(STATUS "AMES macOS target architectures: compiler default (CMake system processor: ${CMAKE_SYSTEM_PROCESSOR})")
        if(CMAKE_SYSTEM_PROCESSOR AND
           NOT CMAKE_SYSTEM_PROCESSOR STREQUAL AMES_PYTHON_ARCH)
            message(WARNING
                "Python is running as ${AMES_PYTHON_ARCH}, but CMake reports "
                "${CMAKE_SYSTEM_PROCESSOR}. Check for mixed ARM64/x86_64 tools. "
                "For native Apple Silicon builds, select an ARM64 Python and "
                "configure a fresh directory with -DCMAKE_OSX_ARCHITECTURES=arm64.")
        endif()
    endif()
endif()

# Prefer normal CMake/config-package discovery. If it fails, try the pybind11
# package installed for the selected interpreter, without overriding an explicit
# pybind11_DIR or pybind11_ROOT hint.
set(_ames_explicit_pybind11_hint FALSE)
if((DEFINED pybind11_DIR AND pybind11_DIR AND
    NOT pybind11_DIR MATCHES "-NOTFOUND$") OR
   (DEFINED pybind11_ROOT AND pybind11_ROOT) OR
   (DEFINED ENV{pybind11_ROOT} AND NOT "$ENV{pybind11_ROOT}" STREQUAL ""))
    set(_ames_explicit_pybind11_hint TRUE)
endif()
find_package(pybind11 CONFIG QUIET)
if(NOT pybind11_FOUND AND NOT _ames_explicit_pybind11_hint)
    execute_process(
        COMMAND "${Python_EXECUTABLE}" -m pybind11 --cmakedir
        OUTPUT_VARIABLE _ames_pybind11_dir
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
        RESULT_VARIABLE _ames_pybind11_status)
    if(_ames_pybind11_status EQUAL 0 AND
       EXISTS "${_ames_pybind11_dir}/pybind11Config.cmake")
        set(pybind11_DIR "${_ames_pybind11_dir}" CACHE PATH
            "pybind11 CMake package for the selected Python" FORCE)
    endif()
endif()
find_package(pybind11 CONFIG REQUIRED)

# Clang on macOS needs an external OpenMP runtime. Preserve the selected Python
# environment's runtime when present; never silently mix its library with headers
# from another prefix. AMES_OPENMP_ROOT is a user hint, not an auto-filled cache.
if(APPLE AND CMAKE_CXX_COMPILER_ID MATCHES "^(AppleClang|Clang)$")
    set(_ames_openmp_root "")
    if(AMES_OPENMP_ROOT)
        set(_ames_openmp_root "${AMES_OPENMP_ROOT}")
    elseif(EXISTS "${AMES_PYTHON_PREFIX}/lib/libomp.dylib")
        set(_ames_openmp_root "${AMES_PYTHON_PREFIX}")
    else()
        find_program(AMES_BREW_EXECUTABLE NAMES brew
            HINTS /opt/homebrew/bin /usr/local/bin)
        if(AMES_BREW_EXECUTABLE)
            execute_process(
                COMMAND "${AMES_BREW_EXECUTABLE}" --prefix libomp
                OUTPUT_VARIABLE _ames_brew_openmp_root
                OUTPUT_STRIP_TRAILING_WHITESPACE
                ERROR_QUIET
                RESULT_VARIABLE _ames_brew_openmp_status)
            if(_ames_brew_openmp_status EQUAL 0 AND _ames_brew_openmp_root)
                set(_ames_openmp_root "${_ames_brew_openmp_root}")
            endif()
        endif()
        if(NOT _ames_openmp_root)
            foreach(_ames_candidate /opt/homebrew/opt/libomp /usr/local/opt/libomp)
                if(EXISTS "${_ames_candidate}/lib/libomp.dylib")
                    set(_ames_openmp_root "${_ames_candidate}")
                    break()
                endif()
            endforeach()
        endif()
    endif()

    if(NOT _ames_openmp_root OR
       NOT EXISTS "${_ames_openmp_root}/lib/libomp.dylib")
        message(FATAL_ERROR
            "OpenMP runtime not found at '${_ames_openmp_root}/lib/libomp.dylib'. "
            "Install the runtime and development headers in the selected Python "
            "environment, or use 'brew install libomp' for a Homebrew setup. "
            "An explicit prefix can be supplied with -DAMES_OPENMP_ROOT=/path/to/prefix.")
    endif()
    if(NOT EXISTS "${_ames_openmp_root}/include/omp.h")
        message(FATAL_ERROR
            "OpenMP library found, but '${_ames_openmp_root}/include/omp.h' is missing. "
            "Install the matching development headers in this same prefix. "
            "AMES will not silently switch to another runtime when headers are missing.")
    endif()

    if(CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
        set(_ames_openmp_flags "-Xpreprocessor -fopenmp")
    else()
        set(_ames_openmp_flags "-fopenmp=libomp")
    endif()
    set(OpenMP_CXX_FLAGS "${_ames_openmp_flags}" CACHE STRING
        "OpenMP flags for the selected macOS Clang compiler" FORCE)
    set(OpenMP_CXX_LIB_NAMES "omp" CACHE STRING "OpenMP runtime names" FORCE)
    set(OpenMP_CXX_INCLUDE_DIR "${_ames_openmp_root}/include" CACHE PATH
        "OpenMP headers for the selected runtime" FORCE)
    set(OpenMP_omp_LIBRARY "${_ames_openmp_root}/lib/libomp.dylib" CACHE FILEPATH
        "AMES OpenMP runtime" FORCE)
    message(STATUS "AMES OpenMP runtime: ${OpenMP_omp_LIBRARY}")
elseif(AMES_OPENMP_ROOT)
    message(WARNING
        "AMES_OPENMP_ROOT is only used with Clang on macOS. "
        "This compiler/platform uses CMake's standard OpenMP discovery.")
endif()

find_package(OpenMP REQUIRED COMPONENTS CXX)
message(STATUS "AMES OpenMP libraries: ${OpenMP_CXX_LIBRARIES}")

# A preset path is not proof of a working OpenMP toolchain. Compile AND link a
# small program using the same imported target that AMES will use. Recheck each
# configure so changed runtime paths cannot reuse a stale successful result.
include(CheckCXXSourceCompiles)
include(CheckCXXCompilerFlag)
include(CMakePushCheckState)
function(ames_check_openmp)
    cmake_push_check_state(RESET)
    set(CMAKE_REQUIRED_LIBRARIES OpenMP::OpenMP_CXX)
    set(CMAKE_TRY_COMPILE_TARGET_TYPE EXECUTABLE)
    unset(AMES_OPENMP_COMPILES_AND_LINKS CACHE)
    check_cxx_source_compiles([=[
        #include <omp.h>
        #ifndef _OPENMP
        #error "OpenMP is not enabled"
        #endif
        int main() {
            int threads = 0;
            #pragma omp parallel reduction(+:threads)
            { threads += 1; }
            return (threads > 0 && omp_get_max_threads() > 0) ? 0 : 1;
        }
    ]=] AMES_OPENMP_COMPILES_AND_LINKS)
    cmake_pop_check_state()
    if(NOT AMES_OPENMP_COMPILES_AND_LINKS)
        message(FATAL_ERROR
            "The OpenMP compile/link check failed. Check the compiler, target "
            "architecture, headers, and runtime library. See CMake's configure "
            "diagnostics in the build directory. On macOS, verify that libomp "
            "contains every architecture requested for AMES.")
    endif()
endfunction()
ames_check_openmp()

pybind11_add_module(AMES
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

target_include_directories(AMES PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_compile_definitions(AMES PRIVATE
    USE_OPENMP
    PHOTOPAIR_USE_OPENMP
    INFO
    DATA_PATH="${CMAKE_CURRENT_SOURCE_DIR}"
)

# pybind11 handles Python-extension linking; do not add an explicit libpython
# link. The imported OpenMP target supplies both compile and link requirements.
target_link_libraries(AMES PRIVATE OpenMP::OpenMP_CXX)

# Preserve the existing symbol-visibility behavior in this build-only revision.
set_target_properties(AMES PROPERTIES
    CXX_VISIBILITY_PRESET default
    VISIBILITY_INLINES_HIDDEN OFF
)

if(AMES_FAST_MATH)
    check_cxx_compiler_flag("-ffast-math" AMES_HAS_FAST_MATH_FLAG)
    if(NOT AMES_HAS_FAST_MATH_FLAG)
        message(FATAL_ERROR
            "This compiler does not support -ffast-math. "
            "Configure with -DAMES_FAST_MATH=OFF or select a supported compiler.")
    endif()
    target_compile_options(AMES PRIVATE -ffast-math)
endif()

if(AMES_NATIVE_ARCH)
    if(CMAKE_CROSSCOMPILING)
        message(FATAL_ERROR "AMES_NATIVE_ARCH must be OFF for cross-compilation.")
    endif()
    if(APPLE)
        list(LENGTH CMAKE_OSX_ARCHITECTURES _ames_arch_count)
        if(_ames_arch_count GREATER 1)
            message(FATAL_ERROR "AMES_NATIVE_ARCH must be OFF for a universal macOS build.")
        endif()
    endif()
    check_cxx_compiler_flag("-march=native" AMES_HAS_NATIVE_ARCH_FLAG)
    if(NOT AMES_HAS_NATIVE_ARCH_FLAG)
        message(FATAL_ERROR
            "This compiler/target does not support -march=native. "
            "Configure with -DAMES_NATIVE_ARCH=OFF.")
    endif()
    target_compile_options(AMES PRIVATE -march=native)
endif()

# Preserve optional externally supplied Boost settings without adding a new
# mandatory dependency or mixing target_link_libraries signatures.
if(Boost_FOUND)
    if(Boost_INCLUDE_DIRS)
        target_include_directories(AMES PRIVATE ${Boost_INCLUDE_DIRS})
    endif()
    if(Boost_LIBRARIES)
        target_link_libraries(AMES PRIVATE ${Boost_LIBRARIES})
    endif()
endif()

message(STATUS "AMES build type: ${CMAKE_BUILD_TYPE} (multi-config generators select it at build time)")
message(STATUS "AMES fast math: ${AMES_FAST_MATH}")
message(STATUS "AMES native CPU optimization: ${AMES_NATIVE_ARCH}")

# This intentionally remains an absolute destination for the selected Python.
install(TARGETS AMES LIBRARY DESTINATION "${Python_SITEARCH}")
