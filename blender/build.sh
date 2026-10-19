#!/bin/bash
set -euxo pipefail

pip download pyscaffolder --dest ./wheels --only-binary=:all: --python-version=3.13 --platform=win_amd64
pip download pyscaffolder --dest ./wheels --only-binary=:all: --python-version=3.13 --platform=manylinux_2_28_x86_64
pip download pyscaffolder --dest ./wheels --only-binary=:all: --python-version=3.13 --platform=macosx_12_0_arm64
pip download pyscaffolder --dest ./wheels --only-binary=:all: --python-version=3.13 --platform=macosx_12_0_x86_64