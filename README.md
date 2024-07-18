# AMES (Astrophysical Multimessenger Emission Simulator)
This code is developed by

Prof. Kohta Murase,
Dr. Bing Theodore Zhang,
et al

## Install
Check cmake and swig:

cmake version: > 3.14

swig version: > 4.0 (brew install swig)
```
mkdir AMES-install
cd AMES-install
cmake /Users/bzhang/Research/Code/GitHub/AMES
make
```
### [Note] Replace above path to your AMES path

### [Note] Please add following two lines in CMakeLists.txt if you are using M1/2 Mac
```
set(CMAKE_OSX_ARCHITECTURES "arm64" CACHE STRING "The list of target architectures to build")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -arch arm64")
```

### [Warnning] There is an issue of install with python virtual environment for Mac OS, e.g., (base)

## Set python file path
If cmake and make is successful, a python file AMES.py will be generated.

Open .zshrc, and add
```
export PYTHONPATH="${PYTHONPATH}:/Users/bzhang/Research/Code/AMES-install"
```
### [Note] Replace above path to your AMES-install path

Then,
```
source ~/.zshrc
```

## Run AMES
```
import AMES
```
### [Warnning] keep the same version of python used in camke.


