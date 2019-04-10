# ChangeLog

## v2.0.0-rc.1

* Wrapping of NiftyReg to allow registration/resampling in SIRF.
* Implemented new `ImageData` hierarchy common to PET and MR. `ImageData` contain geometrical info.
* MR/Gadgetron
  * Added default constructor and set_up to MRAcquisitionModel
  * Implemented sorting of MR images
  * Implemented reading of MR acquisition data from ISMRMRD file
* PET/STIR
  * projectors can now handle subsets (although with a somewhat ugly work-around)
  * added FBP2D, SSRB and the Parallel Level Sets prior
  * added TOF bins dimension to `PETAcquisitionData` (still fixed to have size 1)
* C++ changes
  * Removed `using` statements from the C++ header files
  * Created namespace `sirf`
  * include files are now moved to subdirectories (such as `include/sirf/common`).
  * Modified ObjectHandle type so that it can handle both `std::shared_ptr` and `boost::shared_ptr`.
* Python/MATLAB:
  * `petmr_data_path` is now obsolete. Use `examples_data_path` instead.
* Python:
  * everything is now in a `sirf` module. Use for instance `import sirf.Gadgetron`
* Matlab:
  * in keeping with changes to c++ and python, classes are now called with e.g., `sirf.STIR.obj` instead of `mSTIR.obj`. Aliases can be used to shorten this (e.g., `PET=set_up_PET()` and then `PET.obj`).
* CMake:
  * Updated minimum required version of CMake to 3.9.0.

## v1.1.0

* Various bug fixes and corrections
* `BUILD_STIR_WITH_OPENMP` is now `ON` by default
* Gadgetron data processors check for Gadgetron server crash
* More data files in `SIRF/data/examples/MR`
* Grayscale plotting enabled

## v1.1.0-rc.1

* Created a python `sirf` package (recommended way of importing)
  * aliased `p(Gadgetron|STIR|Utilities) -> sirf.(Gadgetron|STIR|Utilities)`
  * added `setup.py`
  * exposed cmake variable `PYTHON_STRATEGY`. Options:
     * `PYTHONPATH`: prefix `$PYTHONPATH` (default)
     * `SETUP_PY`:   execute `${PYTHON_EXECUTABLE} setup.py install`
     * `CONDA`:      do nothing
* Added `PYTHON_DEST_DIR` variable, which allows the user to select the install destination of the SIRF python modules. `PYTHON_DEST_DIR` is a cached variable which can be updated on the GUI. If `PYTHON_DEST_DIR` is not set, we will install in `${CMAKE_INSTALL_PREFIX}/python`. Likewise for `MATLAB_DEST_DIR`.
* Some improvements to the demos. Note that PET reconstruction demos have somewhat different parameters.
* Implemented PLS Prior
* Implemented 2D Filtered Back Projection

## v1.0.0

* Access to all MR images and acquisition parameters
* All 8 file IO available (PET: Interfile, MR: HDF5)
* PET
  * PETAcquisitionData object creation from scanner name and parameters
  * ListmodeToSinograms converter class, also estimating randoms (from delayed coincidences)
  * Normalization from ECAT8 (Siemens mMR) and attenuation image
  * Build with OpenMP delivers stable and substantially accelerated performance
* More documentation
  * Developer's Guide
  * Doxygen inline documentation (available on CCP PETMR website)
* More tests (now run via CTest), for Python, Matlab and C++.
* Coverage reporting for Python tests done by ctest

## v0.9.2

- fixed version number and avoid confusing with wrong tag v0.9.1

## v0.9.1

- PET data algebra implemented
- Storage scheme (file/memory) management for acquisition data implemented
- Using single precision float Matlab and Python arrays now
- Argument validity checks introduced
- Consistent naming scheme for libraries and modules adopted
- Matlab tests added
- User Guide Appendix on advanced features added
	- Storage scheme management
	- Programming Gadgetron chains
- Specific versions of dependencies (ISMRMRD, Gadgetron, STIR, SIRF) in SuperBuild
- SuperBuild update for Virtual Machine

## v0.9.0

- first release
