# General information

This document currently only gives major releases. Naming of releases in
this document is “feature-based”. However, in practice we will use
“semantic versioning” (only incrementing major release number if there
are backwards compatibility issues).

Detailed (short-term) issues are at

<https://github.com/CCPPETMR/SIRF-SuperBuild/milestones>

<https://github.com/CCPPETMR/SIRF/milestones>

<https://github.com/CCPPETMR/CCPPETMR_VM/milestones>

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

# SIRF 1.x
These updates will be split over a few intermediate releases.

Target date: Q2 2018

  - Software
    
      - Addition of non-TOF scatter
    
      - Partial support for GE Signa PET/MR (PET data only).

  - Improved documentation

  - CMake SuperBuild of SIRF and Windows (Gadgetron not yet on Windows)

  - Small database with phantom data for testing

# SIRF 2.0

Target date: Q4 2018

  - Add major features that didn’t make it into SIRF 1.x

  - Software
    
      - C++ Interface. (enabling a possible move to SWIG for supporting
        other languages)
    
      - Geometric information encoded in Image objects (coregistered PET
        and MR). Reinterpolating to a different grid size, transforming
        using rigid transformations.
    
      - Common Image objects for PET and MR. Therefore the Shape classes
        will work for PET and MR
    
      - PET reconstruction with MR anatomical priors
    
      - MR iterative reconstruction via Gadgetron
    
      - Full support for measured data (Siemens, GE non-TOF). MR only if
        ismrmrd converter available
        
          - MR sequences: TBC
        
          - PET dynamics and gated (separate reconstructions)
    
      - PET TOF support (no scatter)
    
      - Create subsets of acquisition and image data
      
      - Interface to motion estimation software (NiftyReg?) (Flagship)

  - Expanded Testing framework

  - Installers with precompiled software

# SIRF 3.0

Target date: Q2 2019

  - Software:
    
      - Motion-guided reconstruction (Flagship)  
        Spatial only at first, time sync later
    
      - Functions to compute gradients and values of objective functions
    
      - MR reconstruction with PET prior
    
      - Joint PET-MR reconstruction using MATLAB or Python
        tools/toolboxes
    
      - ADMM implementation of joint reconstruction
    
      - Implementation of a few generic optimisation algorithms
    
      - Support for measured data
        
          - MR sequences: (list TBD)
        
          - GE TOF
        
          - Philips?

  - Documentation on developing new functionality and interfaces
    (“engines”)

# Future

  - Code optimization (GPU?)

  - List mode reconstruction (Flagship)

  - non-standard MR sequences

  - Dynamic/gated data with parametric models (Flagship)

  - Joint motion estimation (Flagship)

  - Non-cuboid voxelised images (e.g. blobs, non-Cartesian grids,
    wavelet representations, etc),

  - Integration of other reconstruction packages
