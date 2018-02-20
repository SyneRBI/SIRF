# ChangeLog

## v1.0

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
