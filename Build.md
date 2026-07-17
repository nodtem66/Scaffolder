# Build Reference

## Quick start

```bash
# Configure with default options (binary only)
cmake -S . -B build

# Build the binary
cmake --build build --config Release

# Build and install
cmake --build build --target install
```

## Options (`-D<VAR>=<VALUE>`)

### Build targets

| Option | Default | Description |
|--------|---------|-------------|
| `BUILD_SCAFFOLDER` | `ON` | Build the native `Scaffolder` binary and `Scaffolder.SliceTest` utility. |
| `BUILD_PYSCAFFOLDER` | `OFF` | Build the `PyScaffolder` Python module (requires Python and PyBind11). |
| `BUILD_WITH_OPENMP` | `ON` | Enable OpenMP parallelisation. |

### libigl

| Option | Default | Description |
|--------|---------|-------------|
| `LIBIGL_WITH_OPENGL` | `OFF` | Build libigl with OpenGL support. |
| `LIBIGL_WITH_OPENGL_GLFW` | `OFF` | Build libigl with GLFW window support. |
| `LIBIGL_BUILD_TESTS` | `OFF` | Build libigl unit tests. |
| `LIBIGL_BUILD_TUTORIALS` | `OFF` | Build libigl tutorial examples. |
| `LIBIGL_USE_STATIC_LIBRARY` | `ON` | Link libigl as a static library. |
| `IGL_STATIC_RUNTIME` | `ON` | Use static runtime for libigl. |

### External dependencies

These variables are only consulted when the required external dependency tree is **missing** from the `external/` folder. They control the automatic download fallback. `external.tar.gz` contains snapshot version of libigl, openctm, vcglib, sol2, and pybind11 for building the software. Cmake automatically downloads this tar by default. This provides a conservative release model like Debian which prioritizes long-term stability and reliability over cutting-edge volatility.


| Variable | Default | Description |
|----------|---------|-------------|
| `SCAFFOLDER_EXTERNAL_URL` | `https://github.com/nodtem66/Scaffolder/releases/download/v1.5.3/external-154.tar.gz` | Download URL for the external dependencies archive. |
| `SCAFFOLDER_EXTERNAL_ARCHIVE` | `<root>/scripts/external.tar.gz` | Local path where the archive is saved. |
| `SCAFFOLDER_EXTERNAL_HASH` | `c30df13f296a8a1c421164d7c5c724ce4481fb318d863200190b5ca5c33f53b9` | SHA-256 hash for archive integrity verification. |


### Compiler and platform

| Setting | Scope | Description |
|---------|-------|-------------|
| `CMAKE_CXX_STANDARD` | Project | C++ standard — **17** (hard-set). |
| `CMAKE_OSX_ARCHITECTURES` | macOS | Universal binary architectures — `x86_64;arm64`. |
| `OpenMP_CXX_FLAGS` | macOS + Clang | Preprocessor flags for OpenMP on Apple Clang. |
| `OpenMP_C_FLAGS` | macOS + Clang | Preprocessor flags for OpenMP on Apple Clang. |

## Usage examples

### Build only the CLI binary (default)

```bash
cmake -S . -B build
cmake --build build
```

### Build the Python module only

```bash
cmake -S . -B build_py -DBUILD_SCAFFOLDER=OFF -DBUILD_PYSCAFFOLDER=ON
cmake --build build_py
```

### Disable OpenMP

```bash
cmake -S . -B build -DBUILD_WITH_OPENMP=OFF
cmake --build build
```

### Override the external dependency archive

```bash
cmake -S . -B build \
    -DSCAFFOLDER_EXTERNAL_URL="https://example.com/custom-external.tar.gz" \
    -DSCAFFOLDER_EXTERNAL_HASH="abcdef..."
```

## Build targets

| Target | Description |
|--------|-------------|
| `Scaffolder` | The main CLI binary. |
| `Scaffolder.SliceTest` | Standalone slice-test utility. |
| `PyScaffolder` | Python bindings module (built when `BUILD_PYSCAFFOLDER=ON`). |

## Install components

| Component | Targets | Destination |
|-----------|---------|-------------|
| `scaffolder` | `Scaffolder`, `Scaffolder.SliceTest` | `bin/` |
| `python` | `PyScaffolder` | Standard CMake install prefix (`lib/`, `bin/`, or `Lib/site-packages/` depending on platform). |

## CMake policies

| Policy | Setting | Reason |
|--------|---------|--------|
| `CMP0135` | `NEW` | Allow `FetchContent` to use `URL_HASH` for integrity verification. |
| `CMP0169` | `OLD` | Preserve pre-3.18 `FetchContent` behaviour for `ExternalProject` integration. |

## Required system dependencies

The build will automatically download missing dependencies into `external/` if they are not already present. The following packages are expected:

- **libigl** — geometry processing library
- **sol2** — Lua/C++ binding
- **pybind11** — Python/C++ binding (only for `PyScaffolder`)
- **VCGlib** — mesh processing library
- **OpenCTM** — mesh compression format
- **Eigen3** — linear algebra (bundled with VCGlib)
- **Lua 5.3.5** — scripting language (built from source by the project)
- **Threads** — C++ threading library
- **OpenMP** — shared-memory parallel programming (optional, controlled by `BUILD_WITH_OPENMP`)