#!/bin/bash
set +e
set -u

### Create Makefiles
cmake -G Ninja \
      -DCMAKE_INSTALL_PREFIX=$PREFIX \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_CXX_STANDARD=11 \
      -DBUILD_PYSCAFFOLDER:BOOL=OFF \
      -S . -B build
EXITCODE=$?
if [[ $EXITCODE -ne 0 ]]; then
  [[ -f ${SRC_DIR}/build/CMakeFiles/CMakeOutput.log ]] && cat ${SRC_DIR}/build/CMakeFiles/CMakeOutput.log
  [[ -f ${SRC_DIR}/build/CMakeFiles/CMakeError.log ]] && cat ${SRC_DIR}/build/CMakeFiles/CMakeError.log
  exit $EXITCODE
fi

### Build
cmake --build build -- -j${CPU_COUNT}
EXITCODE=$?
if [[ $EXITCODE -ne 0 ]]; then
  [[ -f ${SRC_DIR}/build/CMakeFiles/CMakeOutput.log ]] && cat ${SRC_DIR}/build/CMakeFiles/CMakeOutput.log
  [[ -f ${SRC_DIR}/build/CMakeFiles/CMakeError.log ]] && cat ${SRC_DIR}/build/CMakeFiles/CMakeError.log
  exit $EXITCODE
fi
set -e

### Install
cmake --install build --component scaffolder

### Test / Check ?
### There is no make check/test