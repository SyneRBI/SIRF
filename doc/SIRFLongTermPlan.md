# General information

This document currently only gives major releases. Naming of releases in
this document is “feature-based”. However, in practice we will use
“semantic versioning” (only incrementing major release number if there
are backwards compatibility issues).

Detailed (short-term) issues are at

- <https://github.com/SyneRBI/SIRF-SuperBuild/milestones>
- <https://github.com/SyneRBI/SIRF/milestones>
- <https://github.com/SyneRBI/SyneRBI_VM/milestones>

Descriptions on past releases are not necessarily complete. Check our
[../CHANGES.md](CHANGES.md) for more information.

# SIRF 0.9
Released 16 Oct 2017.

SIRF 0.9 was our first public release. There will be some interface
changes between 0.9 and 1.0.

# SIRF 1.0 
Released 3 April 2018.

  - Software
      - Common framework for PET and MR reconstruction with MATLAB and
        Python interfaces.
      - Full support for Siemens mMR (both PET and MR).
      - Static data reconstruction, independent for PET and MR
  - Basic testing (including Continuous Integration)
  - User documentation, including inline help
  - Installation options:
      - CMake SuperBuild of SIRF on Linux, MacOS
      - Virtual Machine (VirtualBox) with pre-installed software
      - Docker images with pre-installed software

# SIRF 2.0 
Released 14 May 2019.

  - Software
      - Improved C++ interface with growing similarity with the Python/MATLAB classes.
         - Common `ImageData` objects for PET and MR.
      - Geometric information encoded in `ImageData` objects (coregistered PET
        and MR). Reinterpolating to a different grid size, transforming
        using rigid transformations.
      - Interface to motion estimation software (via [NiftyReg](http://cmictig.cs.ucl.ac.uk/wiki/index.php/NiftyReg)) (Flagship)

      - PET reconstruction with MR anatomical priors.
    - Improved documentation


# SIRF 2.1 
Released 20 Nov 2019.

  - Software
      - Incorporation of the Hybrid Kernel EM method for PET reconstruction with MR info.
      - MR acquisition modelling of 3D Cartesian sequences (with undersampling)
      - Possibility to add DICOM output to an MR reconstruction chain
      - Integration with the [Core Imaging Library (CIL)](https://github.com/vais-ral/CCPi-Framework) (Python-only) for access to general optimisers and regularisation.

# SIRF 2.2
Released May 2020.

  - Software
      - Improvements to resampling, including adjoint operation
      - Ability to interface to SPM registration
      - Improvements in image file format conversion, including LPS information in MR reconstructed images
      - Further integration with the [Core Imaging Library (CIL)](https://github.com/vais-ral/CCPi-Framework) (Python-only) for access to general optimisers and regularisation.

# SIRF 3.0
Released 20 May 2020

  - Drop support for Python 2.
  - Simplify handling of coil images and sensitivities in C++ code.
  - `Registration`: rename `Resample` to `Resampler` and `NiftyResample` to `NiftyResampler`.
  - Restrict MR acquisition data storage scheme to memory only.
  - Add of `sirf.STIR.ScatterEstimation` and `ScatterSimulation` to allow (non-TOF) scatter estimation in PET.
  - Support GE Signa PET/MR reading of listmode data, sinograms, normalisation and randoms.
  - Make [Georg Schramm's parallel (GPU) projector](https://github.com/gschramm/parallelproj) available from SIRF.
  - Implement extraction of the operator representing the linear part of PET acquisition model and computation of its norm.
  - Support registration of 2d images with aladin and f3d.

# SIRF 3.1
Released 24 Jun 2020

  - Support golden-angle radial phase encoding (RPE) if `Gadgetron` toolboxes were found during built.
  - Use `sort_by_time()` for MR acquisitions.
  - Introduce encoding classes that perform the Fourier transformations instead of the `MRAcquisitionModel`.
  - Add constructor for `GadgetronImagesVector` from `MRAcquisitionData`.
  - Methods `range_geometry` and `domain_geometry` of `AcquisitionModel` classes in Python interface, required by CIL algorithms, should obtain data via respective C++ `AcquisitionModel` classes accessors.

# SIRF 3.2
Target: 2022 Q1

  - Replace where possible returning `stir::Succeeded::no` with throwing exception.
  - Remove `__div__` ,  `__idiv__` operators for `DataContainer` classes in Python interface required for Python2.
  - Add `__truediv__` and `__itruediv__` Python3 operators to `DataContainer` classes.
  - Export a CMake config file such that external C++ projects can use SIRF via CMake.
  - During the build step run the executable `ismrmrd_generate_cartesian_shepp_logan` to generate test data compatible with the installed ISMRMRD version.
  - Require ISMRMRD not older than v1.4.2.1.
  - Add conjugation methods to `DataContainer` classes.

# SIRF 3.3
Target: Q4 2022

  - Create subsets of acquisition and image data.
  - C++ functions to compute gradients and values of MR objective function.
  - More MR sequences (list of example sequences that we support TBC).
  - PET TOF support (no scatter; depends on STIR).
  - Implement creation of gated PET sinograms from listmode (depends on STIR).
  - Motion estimation via image registration in SIRF (via connection with CIL).
  - Further additions to C++ Interface coherent with Python and Matlab interface (enabling a possible move to direct use of SWIG).

# SIRF 4.0
Target: Q4 2023

  - Software
      - LPS coordinate system that coincides with the vendor's, including handling of bed position.
      - Full support for measured data (Siemens, GE non-TOF; MR only if ismrmrd converter available).
      - PET dynamics and gated (separate reconstructions), needs support on `DataContainer` classes.
      - Extended SPECT support, including multiple energy windows and scatter (via STIR).
      - Expanded Testing framework.
      - PET list mode reconstruction.
      - Joint motion and reconstruction estimation (via connection with CIL; spatial only at first, time sync later).
      - MR reconstruction with PET prior (via connection with CIL).
      - Joint PET-MR reconstruction using MATLAB or Python tools/toolboxes.
      - Implementation of a few generic optimisation algorithms (to be decided how much of this needs to be in SIRF vs CIL).
  - Sample pipelines for PET and MR reconstruction for static data (based on current scripts)
      - Add error checking of input.
      - Standardise input and output file structure/location.

# Future
  - Installers with precompiled software (conda).
  - Strategy for developing new functionality and interface (“engines”).
  - Additional support for measured data
     - MR sequences: (list TBD)
     - Philips?
  - Example interfaces to Machine Learning framework(s).
  - Further ode optimization.
  - Dynamic/gated data with parametric models.
  - Non-cuboid voxelised images (e.g. blobs, non-Cartesian grids, wavelet representations, etc),
  - Integration of other reconstruction packages.
