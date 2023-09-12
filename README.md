# AMES (Astrophysical Multimessenger Emission Simulator)
This code is developed by

Dr. Kohta Murase,
Dr. Bing Theodore Zhang,
et al

Citations:
Zhang, Murase, Veres & Meszaros, Astrophys.J. 920 (2021) 55

## Build
The current version not support intallation within python virtual environment.

Check cmake and swig:

cmake version: > 3.14

swig version: > 4.0
```
mkdir build && cd build
cmake ..
make
```
If cmake and make is successful, a python file AMES.py will be generated.

## Run
```
import AMES

Follow the examples in the main folder.

```
