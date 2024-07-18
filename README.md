# AMES (Astrophysical Multimessenger Emission Simulator)

**Developed by:**
- Prof. Kohta Murase,
- Dr. Bing Theodore Zhang,
- et al

## Installation

### Prerequisites

Ensure you have the required versions of cmake and swig installed:

- **cmake version:** > 3.14

- **swig version:** > 4.0 (install via `brew install swig')

### Installation Steps

1. Create and navigate to a directory for the AMES installation:
   ```
   mkdir AMES-install
   cd AMES-install
   ```

2. Run cmake and make:
   ```
   cmake /Users/bzhang/Research/Code/GitHub/AMES
   make
   ```

#### Note: Replace `/Users/bzhang/Research/Code/GitHub/AMES` with the actual path to your AMES directory.

### Additional Steps for M1/M2 Mac Users
If you are using an M1 or M2 Mac, add the following two lines to `CMakeLists.txt`:

```
set(CMAKE_OSX_ARCHITECTURES "arm64" CACHE STRING "The list of target architectures to build")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -arch arm64")
```

### Warning for Mac OS Users
There is an issue with installing AMES within a Python virtual environment on Mac OS (e.g., (base)). It is recommended to use the system Python environment.

## Set python file path
If cmake and make are successful, a Python file `AMES.py` will be generated.

1, Open `.zshrc` file and add the following line:
```
export PYTHONPATH="${PYTHONPATH}:/Users/bzhang/Research/Code/AMES-install"
```

#### Note: Replace `/Users/bzhang/Research/Code/AMES-install` with the actual path to your AMES-install directory.

2, Reload your `.zshrc`:
```
source ~/.zshrc
```

## Running AMES

To use AMES in your Python environment, simply import it:
```
import AMES
```
### Warnning: 

Ensure that you use the same Python version that was used during the cmake process.

