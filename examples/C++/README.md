This folder contains example C++ files that use SIRF, together with
associated CMake files.

SIRF exports a `SIRFConfig.cmake` etc. This defines targets for every
SIRF library, e.g. `sirf::common`, `sirf::Reg` etc
We can therefore use
```cmake
find_package(SIRF)
```
When running CMake, you might have to set SIRF_DIR to let CMake find the exported
`SIRFConfig`.


To build this project, follow the usual CMake steps. For instance, when using
the command version of CMake:
```sh
mkdir build
cd build
cmake /wherever/is/the/source/examples/C++ \
   -DSIRF_DIR:PATH=/wherever/is/the/SIRF/install/lib/cmake/SIRF-3.1
```
