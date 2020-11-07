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

SIRF 0.9 was our first public release. There will be some interface
changes between 0.9 and 1.0.

# SIRF 1.0 (released 3 April 2018)

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

# SIRF 2.0 (released 14 May 2019)

  - Software
      - Improved C++ interface with growing similarity with the Python/MATLAB classes.
         - Common `ImageData` objects for PET and MR.
      - Geometric information encoded in `ImageData` objects (coregistered PET
        and MR). Reinterpolating to a different grid size, transforming
        using rigid transformations.
      - Interface to motion estimation software (via [NiftyReg](http://cmictig.cs.ucl.ac.uk/wiki/index.php/NiftyReg)) (Flagship)

      - PET reconstruction with MR anatomical priors.
    - Improved documentation


# SIRF 2.1 (released 20 Nov 2019)

  - Software
      - Incorporation of the Hybrid Kernel EM method for PET reconstruction with MR info.
      - MR acquisition modelling of 3D Cartesian sequences (with undersampling)
      - Possibility to add DICOM output to an MR reconstruction chain
      - Integration with the [Core Imaging Library (CIL)](https://github.com/vais-ral/CCPi-Framework) (Python-only) for access to general optimisers and regularisation.

# SIRF 2.2 (released May 2020)

  - Software
      - Improvements to resampling, including adjoint operation
      - Ability to interface to SPM registration
      - Improvements in image file format conversion, including LPS information in MR reconstructed images
      - Further integration with the [Core Imaging Library (CIL)](https://github.com/vais-ral/CCPi-Framework) (Python-only) for access to general optimisers and regularisation.

# SIRF 2.3
Target date : Q4 2020

  - Software
      - Addition of non-TOF scatter
      - Partial support for GE Signa PET/MR (PET data only).
      - Extension of HKEM to use multiple anatomical images as prior
      - MR iterative reconstruction via connection with CIL 
            - improvements on acquisition model (non-Cartesian sequences)
            - selecting of MR aquisition data
      - Motion Compensate Image Reconstruction (MCIR) via connection with CIL, both PET and MR
      - Initial SPECT support (via STIR)
 - optional: CMake SuperBuild of SIRF on Windows (Gadgetron not yet on Windows)
      
# SIRF 3.0
Target date: Q2 2021

  - Add major features that didn’t make it into SIRF 2.x
  - Software
      - LPS coordinate system that coincides with the vendor's, including handling of bed position.
      - Further additions to C++ Interface. (enabling a possible move to SWIG for supporting
        other languages)
      - Functions to compute gradients and values of MR objective function (independent of CIL?).
      - Full support for measured data (Siemens, GE non-TOF). MR only if
        ismrmrd converter available
      - MR sequences: (list of example sequences that we support, TBC (includes standard 2D radial sequences))
      - creation of gated PET sinograms from listmode
      - PET dynamics and gated (separate reconstructions), needs support on `DataContainers`
      - PET TOF support (no scatter) (depends on STIR)
      - Extended SPECT support, including multiple energy windows and scatter (via STIR)
      - Create subsets of acquisition and image data
      - Motion estimation via image registration in SIRF (via connection with CIL)
      - Expanded Testing framework
  - Sample pipelines for PET and MR reconstruction for static data (based on current scripts)
      - Add error checking of input
      - Standardise input and output file structure/location
  - Installers with precompiled software (conda)
  - Strategy for developing new functionality and interfaces
    (“engines”)
    
  Proposed to drop:
      - MR iterative reconstruction via Gadgetron gadgets (e.g. iterative SENSE). Reason: superseded by CIL functionality and probably a lot of work

# SIRF 4.0
Target date: Q4 2021

  - Software:
      - PET List mode reconstruction
      - Joint motion and reconstruction estimation (via connection with CIL)
        Spatial only at first, time sync later
      - MR reconstruction with PET prior (via connection with CIL)
      - Joint PET-MR reconstruction using MATLAB or Python
        tools/toolboxes
      - Implementation of a few generic optimisation algorithms (to be decided how much of this needs to be in SIRF vs CIL)
      - Additional support for measured data
         - MR sequences: (list TBD)
         - Philips?


# Future
  - Code optimization (GPU?)
  - Dynamic/gated data with parametric models
  - Non-cuboid voxelised images (e.g. blobs, non-Cartesian grids,
    wavelet representations, etc),
  - Integration of other reconstruction packages
