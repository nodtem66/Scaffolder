@echo on

cmake %CMAKE_ARGS% -GNinja ^
  -DCMAKE_PREFIX_PATH=%LIBRARY_PREFIX% ^
  -DCMAKE_INSTALL_PREFIX=%LIBRARY_PREFIX% ^
  -DCMAKE_BUILD_TYPE=Release ^
  -DSCAFFOLDER_BUILD_PYTHON=OFF ^
  -DVERSION=%SCAFFOLDER_VERSION% ^
  -S %SRC_DIR% ^
  -B build

cmake --build build --parallel %CPU_COUNT%
cmake --install build