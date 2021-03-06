#!/bin/bash
set -eu

echo -n "Create external.tar.gz    "
if [ ! -e external.tar.gz ]; then
  tar --directory=../ \
      --exclude=external/sol2/docs \
      --exclude=external/examples/docs \
      --exclude=external/tests/scripts \
      --exclude=external/tests/docs \
      --exclude=external/tbb/doc \
      --exclude=external/tbb/examples \
      --exclude=external/tbb/test \
      --exclude=external/vcglib/docs \
      --exclude=external/vcglib/apps \
      --exclude=external/libigl/tutorial \
      --exclude=external/libigl/tests \
      --exclude=external/libigl/external/.cache \
      --exclude=external/pybind11/docs \
      --exclude=external/OpenCTM-1.0.3/tools \
      --exclude=external/OpenCTM-1.0.3/plugins \
      --exclude=external/OpenCTM-1.0.3/bindings \
      --exclude-vcs --exclude-vcs-ignores \
      -czf external.tar.gz external && echo "[OK]"
else
  echo "[Exists]"
fi