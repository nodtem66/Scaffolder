setlocal EnableDelayedExpansion

rem Create Makefiles
cmake -G Ninja ^
      -DCMAKE_INSTALL_PREFIX:PATH="%PREFIX%" ^
      -DCMAKE_BUILD_TYPE:STRING=Release ^
      -DCMAKE_CXX_STANDARD=11 ^
      -S . -B build
if errorlevel 1 goto DEBUG

rem Build
cmake --build build
if errorlevel 1 goto DEBUG

rem Install
cmake --install build --component scaffolder 
if errorlevel 1 goto DEBUG
exit 0

:DEBUG
type build\\CMakeFiles\\CMakeOutput.log
type build\\CMakeFiles\\CMakeError.log
exit 1