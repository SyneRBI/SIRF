# ChangeLog

## v3.0.0
### Backwards incompatible changes
* STIR version 4.1.0 is now required.
* Python 2 is no longer supported. Most code might still work, but we do not check. A warning is written when the Python version found is 2. This will be changed to `FATAL_ERROR` at a later stage. 
* Handling of coil images and sensitivities in C++ code simplified by inheriting CoilImagesVector from GadgetronImagesVector and replacing CoilSensitivitiesAsImages with CoilSensitivitiesVector, also inheriting from GadgetronImagesVector. All methods of CoilImagesVector and CoilSensitivitiesVector other than those inherited from GadgetronImagesVector are no longer supported except methods named compute(), which are renamed to calculate().

### Deprecations (will be errors in SIRF 4.0)
* `Registration`: renamed `Resample` to `Resampler` and `NiftyResample` to `NiftyResampler`. Old names are now deprecated but should still work.
* STIR `AcquisitionModel` `forward`, `direct`, `backward` and `adjoint` signatures have changed in Python. Subset information should now be set via `num_subsets` and `subset_num` members. The `forward` and `backward` members can still be called with the previous syntax but this will be removed in a later version.
Note that default values of `num_subsets` and `subset_num` are 0 and 1 respectively, such that default behaviour is default behaviour (i.e. process all data) is unchanged.
* MR acquisition data storage scheme restricted to memory only (a message will be printed but no error thrown)
* Use CMake variable names from `find_package(Python)` which are available with CMake 3.12+. SIRF CMake files will accept both `Python_EXECUTABLE` or `PYTHON_EXECUTABLE`, for the latter it will send a deprecation warning.

### New features
* PET
  - Addition of `sirf.STIR.ScatterEstimation` and `ScatterSimulation` to allow (non-TOF) scatter estimation in PET
  - GE Signa PET/MR reading of listmode data, sinograms, normalisation and randoms support added.
  - If STIR is at least version 5 or built from the master branch, [Georg Schramm's parallel (computing) projector](https://github.com/gschramm/parallelproj proj) is now made available from SIRF (use `AcquisitionModelUsingParallelproj`). This uses Joseph interpolation, but importantly can use your GPU (if CUDA was found during building).
  - Implemented extraction of the operator representing the linear part of PET acquisition model and computation of its norm.
  - When adding a shape to a `sirf.STIR.ImageData`, optionally give the number of times to sample a voxel. This is useful when the shape partially - but not completely - fills a voxel.
  - If `storage_scheme` is set to `memory`, `PETAcquisitionData` allows direct modification, whereas before a copy would need to be created first. (Internally, it uses STIR `ProjDataInMemory`, instead of `ProjDataFromStream`).
* Registration
  - Registration of 2d images is now supported with aladin and f3d. 
* examples data:
  - Installs `examples`, `data` and `doc` to the install directory, i.e. `${CMAKE_INSTALL_PREFIX}/share/SIRF-<version_major>.<version_minor>` directory.
  - If the `SIRF_DATA_PATH` environment variable is set, `examples_data_path` will search for the examples data there, or in `SIRF_INSTALL_PATH/share/SIRF-<version_major>.<version_minor>/data` directory. In MATLAB, the `example_data_path` function has the version set by CMake at install time.
* Other Python features:
  - Define `__version__` in `sirf` python package.
  - Added implementation of division and multiplication for `NiftiImageData`.
  - Data validity checks return `NotImplemented` instead of throwing error, opening the door for future implementations of operations on data.

### Other changes
* When registering, internally the forward displacement is no longer stored, replaced by the forward deformation. The inverse is no longer stored, and is calculated as needed.
* `PETAcquisitionData.axpby` now uses STIR's `axpby` and is therefore faster.
* Speed-up in `stir::AcquisitionDataInMemory` of `as_array`, `fill`, `dot`, `norm`, etc. (by using STIR iterators).
* Added common Python `DataContainer` algebra unit tests for all `DataContainer` inherited classes.
* Continuous Integration now uses Github Actions. Travis-CI has been dropped.
* New `CMake` option `BUILD_DOCUMENTATION` to use doxygen to build C++ documentation.
It will be installed in the `share/SIRF-version/doc/doxygen`.

### Bug fixes
* Python `fill` method in MR `DataContainer` accepts `numpy` array, number or `DataContainer`.
* `get_index_to_physical_point_matrix()` returned a wrong matrix in MATLAB and Python.
* path manipulation of `examples_data_path` now should work for any platform, not just linux.

## v2.2.0
* Changed CCP PETMR to SyneRBI
* updates to steepest ascent demo
* STIR.AcquisitionData.get_info() returns a string that describes the scanner etc
* documentation fixes/additions

## v2.2.0-rc.1

* A passthrough for both the maximum and minimum relative change during OSMAPOSL reconstruction has been added.
* We have now corrected the geometrical information of `.h5` images (coming from ISMRMRD and Gadgetron). This means we can now convert them to other SIRF image types (e.g., `NiftiImageData` and `STIRImageData`). This is necessary for any kind of synergistic reconstruction. Further, to the best of our knowledge, this is the first ISMRMRD to NIfTI converter out there!
* The adjoint transformation has now been implemented for `NiftyResample` through the wrapping of NiftyMoMo.
* The following methods have been added to C++, python and matlab NiftyResample:
	* `out = forward(in)`
	* `forward(out, in)`
	* `out = adjoint(in)`
	* `adjoint(out, in)`
	* `out = backward(in)` <- alias for adjoint
	* `backward(out, in)` <- alias for adjoint
* Inverse deformation images. Inverse displacements are also possible by converting to and from deformations.
* NiftyPET projector wrapped (if STIR is built with NiftyPET)
* Added `set_image_data_processor` to `PETAcquisitionModel`.  This allows for instance image-based PSF modelling.
* Resampling of complex images.
* SPM registration wrapping (only SPM12 tested). If `Matlab` and `SPM` are present, the SPM wrapper is available from `C++`, `Matlab` and `Python`.
* Support for registering multiple floating images has been added. This is only available for certain algorithms (currently only `SPM`). There are therefore new methods `add_floating_image` and `clear_floating_images` on top of the original `set_floating_image`. Methods extracting the results of registrations can now be called with an index (`get_output(idx = 0)`, `get_transformation_matrix_forward(idx = 0)`, etc.). This index defaults to the first to maintain backwards compatibility.
* Ability to pad `NiftiImageData`, e.g., `a.pad([10,10,0],[10,10,0])` to add 10 voxels to the minimum and maximum of the x- and y-directions.
* Ability to set and get STIR verbosity from python.
* Save STIR images using a parameter file (e.g., for saving as `.nii`)
* Default F3d to using non-symmetric version (previously, symmetric was used). Option to use the symmetric in C++, but currently exposed to python and matlab as we suspect there is an upstream bug there.

## v2.1.0

* PET/STIR
	* Interfaced HKEM into SIRF
	* Interfaced SeparableGaussianImageFilter into SIRF
* MR/Gadgetron
	* Added DICOM-writing gadgets for MR images output
	* Added few Gadgetron GPU gadgets to SIRF gadget library
	* Enabled handling of 3D slices of MR images by switching to 3D FFT
* Python
	* Switched to new class style
	* Introduced contiguity checks of filled data
* CIL/SIRF Compatibility
     * added methods to AcquisitionData, ImageData and AcquisitionModel to be compatible with
       CCPi's Core Imaging Library (CIL)

## v2.0.0

* Set CMake policy CMP0079.
* Use `swig_add_library` instead of `swig_add_module`.
* Averaging of rigid transformation matrices via quaternions (and therefore a quaternion class).
* Arrays of SIRF objects can be passed from the Python and Matlab interfaces to the C++ level (e.g., averaging a large number of matrices) via the DataHandleVector class. This is an internal class that should not be used. Simply pass a native array of objects and SIRF will convert to the DataHandleVector class if necessary.
* Image data role checks in MRAcquisitionModel introduced.
* Corrected ISMRMRD acquisition sorting.
* Added PhysioInterpolationGadget and FatWaterGadget to SIRF gadgets library.
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
  * Doxygen inline documentation (available on SyneRBI website)
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
