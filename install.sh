rm -rf CMakeFiles
rm CMakeCache.txt
rm *AMES*
rm Makefile
export AMES_DIR="$PWD"
CMAKE_PREFIX_PATH=$AMES_DIR cmake -DCMAKE_INSTALL_PREFIX=$AMES_DIR ../GitHub/AMES
make
make install
