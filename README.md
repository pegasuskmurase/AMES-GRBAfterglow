# AMES (Astrophysical Multimessenger Emission Simulator)
This code is developed by

Dr. Kohta Murase, Dr. Bing Theodore Zhang, et al

Citations: Zhang, Murase, Veres & Meszaros, Astrophys.J. 920 (2021) 55

## Install
Check cmake and swig:

cmake version: > 3.14

swig version: > 4.0 (brew install swig)
```
mkdir AMES-install
cmake /Users/bzhang/Research/Code/GitHub/AMES
make
```
### [Note] Replace above path to your AMES path

### [Warnning] Please install without python virtual environment, e.g., (base)

## Set python file path
If cmake and make is successful, a python file AMES.py will be generated.

Open .zshrc, and add
```
export PYTHONPATH="${PYTHONPATH}:/Users/bzhang/Research/Code/AMES-install"
```
Then,
```
source ~/.zshrc
```

## Run AMES
```
import AMES
```
### [Warnning] keep the same version of python used in camke.

```
