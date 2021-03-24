# Installation instructions
To compile and install `hwcart` you need `CMake`, a working MPI environment,
and an installation of the `hwloc` library and headers.

First, clone the repository:

```
git clone git@github.com:NordicHPC/hwcart.git
```

Create the build directory

```
mkdir hwcart/build
cd hwcart/build
```

Configure using `CMake`:
```
cmake ..
```

You can set the following `CMake` variables:

```
HWCART_USE_FORTRAN	(default OFF) Choose whether to build Fortran interface
HWCART_USE_HWLOC	(default ON)  Choose whether to use hwloc (OpenMPI is required if hwloc is not used)
HWCART_USE_TESTS	(default OFF) Choose whether to build tests and examples
```
Note that if the OpenMPI environment is used, it is possible to compile and use `hwcart` without `hwloc`.
In that case `hwcart` will use the OpenMPI hardware locality features instead.

To compile and install:
```
make install
```
